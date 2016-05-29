

classdef exactMergeSolver < handle
    % The class for the exact convex merge solver. 
    % - It first uses a convex program to compute the exact internal
    %   boundary flows at the junction, hence decoupling the PDEs. 
    % - Then it exploits the lax-hopf formula for computing the exact
    %   solutions at any time.
    % We refer to the associated paper for the details.
    % Yanning Li, May 26, 2016
    % network:
    % ---- hwy: a struct, containing information for each hwy link
    % ---- junc: a struct, containing the informatino for each junc, including
    %            the time grid at the junction
    % ---- initial and boundary conditions
    % ---- other auxiliary properties
    
    properties (Access = public)
        
        network_junc;   % struct, .(juncStr) contains junction info
        network_hwy;    % struct, .(linkStr), link info: .IC, .X_grid, .X_grid_cum, .para_~
                        
        junc_labels;    % array, list of junction integer labels
        link_labels;    % array, list of link integer labels
        num_juncs;      % number of junctions in the network
        num_links;      % number of links in the network
        
        t_start;        % the start time
        t_end;          % the end time
        
        M;              % a struct for saving the vehicle ID solution
        k;              % a struct for saving the nonlinearly normalized density in each cell
        
        
        % The following properties are for the optimization program.
        Aineq;      % a matrix with all inequality constraints, Aineq * x <= bineq
        bineq;      % a column vector, Aineq * x <= bineq
        lb;     % column vector, lower bound of decision variables
        ub;     % column vector, upper bound of decision variables
        ctype;  % column vector, type of decision variables
        f;      % column vector, minimize f*x
        dv_index;     % the decision variable index. .(linkStr)=[start,end]
        dv_index_max; % the total number of decision variables
        x;      % x is the solution
        
        t_ref;  % the refernce point for checking the exact solutions
        
        max_search_depth;   % the maximum search depth
        
    end
    
    methods (Access = public)
        %===============================================================
        function self = exactMergeSolver(t_start, t_end)
            % set time horizon
            self.t_start = t_start;
            self.t_end = t_end;
            
            % Construct an empty net object
            self.network_junc = struct;
            self.network_hwy = struct;
            
            self.junc_labels = [];
            self.link_labels = [];
            self.num_links = 0;
            self.num_juncs = 0;
            
            % Here is the state on each link, 
            % self.M.linkstr = [], with dim num_cell+1 x num_step+1
            % discretized iniital condition is saved in the first column
            % left top corner is (0,0) = 0
            self.M = struct;
            
            % Here is the density.
            % Just for visualization. Do not count in the computation time.
            self.k = struct;
            
            % For solving the convex program
            self.Aineq = [];
            self.bineq = [];
            self.lb = zeros(0,1);
            self.ub = zeros(0,1);
            self.ctype = '';
            
            self.f = [];
            self.dv_index = struct;
            
            % default maximum serach depth is 50
            self.max_search_depth = 50;
            
        end
        
        
        %===============================================================
        function addJunc(self, junc, inlabel, outlabel, type_junc, ratio)
            % Create a junction. Only support connection, merge, and diverge
            % input:
            %       junc: int, assign a junction id
            %       inlabel/outlabel: numerical vectors containing the labels of links
            %       type_junc: strings 'connection', 'merge', 'diverge',
            %           'onrampjunc', 'offrampjunc'. The difference between
            %           'onrampjunc' and 'merge' is that 'onrampjunc' has one
            %           'onramp' which is controllable. Same for 'offrampjunc'
            %       ratio: 2x1 numerical priority or distribution coefficients.
            %           e.g. [0.4; 0.6]
            
            % make sure no duplicated label
            if any(self.junc_labels == junc)
                error('ERROR: Duplicated junc label: %d.\n', junc)
            end
            
            juncStr = sprintf('junc_%d', junc);
            self.network_junc.(juncStr).inlabel = inlabel;
            self.network_junc.(juncStr).outlabel = outlabel;
            self.network_junc.(juncStr).type_junc = type_junc;
            self.network_junc.(juncStr).ratio = ratio;
            
            % initialize the time duration and the time grid
            self.network_junc.(juncStr).T = [];
            self.network_junc.(juncStr).T_cum = [];
            
            self.junc_labels = [self.junc_labels; junc];
            self.num_juncs = self.num_juncs + 1;
        end
        
        
        %===============================================================
        function addLink(self, link, para, num_lanes, lengthKM, linkType)
            % Add a link to the network.
            % Link is refered as a segment on the highway. The highway should
            % be divided to links at locations where sensors are deployed, or
            % the property of the road chagnes (drops of lanes).
            % input:
            %       link: numerical label
            %       para: a struct containing road properties
            %       num_lanes: int, the number of lanes
            %       lengthKM: float, length of the link in km
            %       linkType: string, 'freeway', 'onramp', 'offramp'
            
            % make sure no duplicated label
            if any(self.link_labels == link)
                error('ERROR: Duplicated link label.\n')
            end
            
            linkStr = sprintf('link_%d',link);
            
            % set the road parameters
            self.network_hwy.(linkStr).para_vf = para.vf;
            self.network_hwy.(linkStr).para_w = para.w;
            self.network_hwy.(linkStr).para_kc = num_lanes*para.kc_pl; %kc per lane
            self.network_hwy.(linkStr).para_km = num_lanes*para.km_pl;
            self.network_hwy.(linkStr).para_qmax = num_lanes*para.qmax_pl;
            self.network_hwy.(linkStr).para_vmin = para.v_min;
            self.network_hwy.(linkStr).para_vmax = para.v_max;
            
            self.network_hwy.(linkStr).para_postkm = lengthKM;
            self.network_hwy.(linkStr).para_linktype = linkType;
            
            % save the label
            self.link_labels = [self.link_labels; link];
            self.num_links = self.num_links + 1;
            
        end
        
        
        %===============================================================
        function setInitialConForLink(self, link, rho_ini)
            % Set initial condition of each link
            % input:
            %       link: int, link label
            %       rho_ini: struct,
            %               .(IC) num_segments x 1 float, veh/m
            %               .(X_grid_cum) num_segments+1 x 1 float, in m
            
            if ~any(self.link_labels == link)
                error('ERROR: Link %d not defined in the network\n', link);
            end
            
            num_seg = length(rho_ini.IC); 
            linkStr = sprintf('link_%d',link);
            
            % if X_grid_cum is empty, then even discretization
            if ~isfield(rho_ini, 'X_grid_cum') || isempty(rho_ini.X_grid_cum)
                
                tmp_dx = self.network_hwy.(linkStr).para_postkm*1000/num_seg;
                self.network_hwy.(linkStr).X_grid_cum = (0:num_seg)'*tmp_dx;
                
            else
                self.network_hwy.(linkStr).X_grid_cum = rho_ini.X_grid_cum;
            end
            
            self.network_hwy.(linkStr).IC = self.columnize(rho_ini.IC);

        end
        
        
        %===============================================================
        function setInitialCon(self, init_condition)
            % Set initial condition for all links
            % input:
            %       link: int, link label
            %       rho_ini: struct,
            %               .(linkStr).(IC) num_segments x 1 float
            %               .(linkStr).(X_grid_cum) num_segments+1 x 1 float, in m
            
            % if is nan, then evenly discretize each link into cells with
            % length around 200 m. Data will be set as NaN.
            if isempty(init_condition)
                
                for link = self.link_labels'
                    
                    linkStr = sprintf('link_%d', link);
                    % the number of segments of around 200 m length
                    num_seg = ceil(self.network_hwy.(linkStr).para_postkm*1000/200);
                    
                    tmp_dx = self.network_hwy.(linkStr).para_postkm*1000/num_seg;
                    self.network_hwy.(linkStr).X_grid_cum = (0:num_seg)'*tmp_dx;
                    
                    % set the density data as NaN
                    self.network_hwy.(linkStr).IC = ones(num_seg,1)*NaN;
                end
                return
            end
            
            % Otherwise, update the initial conditions
            fields = fieldnames(init_condition);
            
            % warning: the initial condition of some links is not provided
            if length(fields) ~= self.num_links
                warning('WARNING: the initial condition is not completely updated.\n')
            end
            
            % update each link
            for i = 1:length(fields)
                
                linkStr = fields{i};
                
                self.network_hwy.(linkStr).IC =...
                    init_condition.(linkStr).IC;
                num_seg = length(init_condition.(linkStr).IC); 
                % if grid provided, then set;
                % otherwise, idiscretize the space evenly.
                if ~isfield(init_condition.(linkStr), 'X_grid_cum') ||...
                        isempty(init_condition.(linkStr).X_grid_cum)
                
                    tmp_dx = self.network_hwy.(linkStr).para_postkm*1000/num_seg;
                    self.network_hwy.(linkStr).X_grid_cum = (0:num_seg)'*tmp_dx;
                
                else
                    self.network_hwy.(linkStr).X_grid_cum = ...
                        init_condition.(linkStr).X_grid_cum;
                end
                
            end
            
            
            
        end

        
        %===============================================================
        function setBoundaryConForLink(self, link, q_in, q_out, T_in, T_out)
            % Set boundary condition of each link
            % input:
            %       link: int, the link label
            %       q_in: num_in_steps x 1 float
            %       q_out: num_out_steps x 1 float
            %       T_in: num_in_steps x 1 float durations of each step
            %       T_out: num_out_steps x 1 float durations of each step
            % If q_in or q_out is [], then it is decision variable to be
            % copmuted. All values will be set as NaN. The corresponding 
            % T_in, T_out should also be [].
            
            if ~any(self.link_labels == link)
                error('ERROR: Link %d not defined in the network.\n', link);
            end
            
            % the flow data and time must have the same length
            if (~isempty(q_in) && length(q_in)~=length(T_in)) ||...
                    (~isempty(q_out) && length(q_out)~=length(T_out))
                error('Data length must be equal to its time discretization');
            end
            
            linkStr = sprintf('link_%d',link);
            % set the upstream boundary condition if not empty
            if ~isempty(q_in)
                self.network_hwy.(linkStr).T_us = T_in;
                self.network_hwy.(linkStr).T_us_cum = [0; cumsum(T_in)];
                self.network_hwy.(linkStr).BC_us = self.columnize(q_in);
            else
                self.network_hwy.(linkStr).T_us = [];
                self.network_hwy.(linkStr).T_us_cum = [];
                self.network_hwy.(linkStr).BC_us = NaN;
            end
            
            if ~isempty(q_out)
                self.network_hwy.(linkStr).T_ds = T_out;
                self.network_hwy.(linkStr).T_ds_cum = [0; cumsum(T_out)];
                self.network_hwy.(linkStr).BC_ds = self.columnize(q_out);
            else
                self.network_hwy.(linkStr).T_ds = [];
                self.network_hwy.(linkStr).T_ds_cum = [];
                self.network_hwy.(linkStr).BC_ds = NaN;
            end
            
           
        end
        
        
        %===============================================================
        function setBoundaryCon(self, boundary_data)
            % Set boundary condition for all links
            % input:
            %       boundary_data: struct
            %           .(linkStr).T_us, T_us_cum, q_us, v_us; same for ds
            
            % Use the setBoundaryConForLink for each link
            for link  = self.link_labels'
                
                linkStr = sprintf('link_%d',link);
                
                if isfield( boundary_data, linkStr )
                    self.setBoundaryConForLink( link, ...
                                            boundary_data.(linkStr).BC_us,... 
                                            boundary_data.(linkStr).BC_ds,...
                                            boundary_data.(linkStr).T_us,...
                                            boundary_data.(linkStr).T_ds);
                else
                    error('boundary data not complete for all links\n')
                end
            end
           
        end
        
        
        %===============================================================
        function computeExactInternalBoundaryFlows(self, T_junc_init, tol)
            % This function computes the exact internal boudnary flows at
            % the junction using a convex program and the bi-section
            % algorithm
            % input:
            %       T_junc_init, the initial grid at the junction.
            %       tol, number of vehicles, the tolerance for a solution
            %           to be considered as the exact solution
            
            % -------------------------------------------------------
            % Only one merge junction is supported in this version
            if self.num_juncs~=1
                error('exactMergeSolver only supports one merge junction on the network.')
            end
            juncStr = sprintf('junc_%d', self.junc_labels);
            if ~strcmp(self.network_junc.(juncStr).type_junc, 'merge')
                error('exactMergeSolver only supports one merge junction on the network.')
            end
            
            % set the initial grid
            T_junc = T_junc_init;
            
            % -------------------------------------------------------
            % iteratively update the junction grid to get the exact
            % solution
            gotExact = false;
            loopCounter = 0;
            while gotExact == false && loopCounter <=50
                
                loopCounter = loopCounter + 1;
                
                % update the grid
                self.setGridAtJunc(self.junc_labels, T_junc);
                
                % set constratins
                self.setConstraintsAtMerge();
                
                % apply entropy condition at the junction
                self.applyEntropyConAtMerge('all');
                
                % compute the solution
                [self.x, fval, exitflag, output] = self.solveProgram();
                
                [gotExact, steps] = self.checkSolution(tol);

                if gotExact == false
                    T_all = self.updateTimeGrid(steps);
                    T_junc = T_all.(juncStr).T;
                end
                
            end
            
            fprintf('\nGot exact solutions at the junction.\n')
            
            % -------------------------------------------------------
            % Now set the exact solution and associated boundary grid to
            % the boundary conditions of the links
            thru = [];
            for link = self.network_junc.(juncStr).inlabel'
                linkStr = sprintf('link_%d', link);
                self.network_hwy.(linkStr).BC_ds = self.x(...
                    self.dv_index.(linkStr)(1):self.dv_index.(linkStr)(2) );
                if isempty(thru)
                    thru = self.network_hwy.(linkStr).BC_ds;
                else
                    thru = thru + self.network_hwy.(linkStr).BC_ds;
                end
            end
            % set the downstream link inflow as the throughput
            linkStr = sprintf('link_%d', self.network_junc.(juncStr).outlabel);
            self.network_hwy.(linkStr).BC_us = thru;
            
            
        end
        
        
        %===============================================================
        function computeAllSolutionsOnGrid(self, dt, dx)
            % This function computes the solution M(t,x) at all the grid
            % points with resolution dt, dx.
            % input:
            %       dt, dx; the grid resolution, second, meter
            % output:
            %       computed solution will be saved in self.M.(linkStr)
            
            % fundamental diagram for calling single link HJ PDE solver
            fd = struct;
            
            for link = self.link_labels'
                
                linkStr = sprintf('link_%d',link);
                
                % extract the computed result 
                q_in = self.network_hwy.(linkStr).BC_us;
                q_out = self.network_hwy.(linkStr).BC_ds;
                p_ini = self.network_hwy.(linkStr).IC;
                
                % fundamental diagram
                fd.(linkStr) = LH_Tfd(self.network_hwy.(linkStr).para_vf,...
                    self.network_hwy.(linkStr).para_w,...
                    self.network_hwy.(linkStr).para_km);
                
                us_position = 0;
                ds_position = self.network_hwy.(linkStr).para_postkm*1000;
                pbEnv = LH_general(fd.(linkStr),us_position,...
                                    ds_position);
                
                %===========================================
                % extract initial segment vector
                ini_seg = self.network_hwy.(linkStr).X_grid_cum;
                
                % Berkeley toolbox need ini_seg and p_ini to be a row
                % vector
                if ~isrow(ini_seg)
                    ini_seg = ini_seg';
                end
                if ~isrow(p_ini)
                    p_ini = p_ini';
                end
                
                pbEnv.setIniDens(ini_seg,p_ini);
                
                % set upstream boundary condition
                time_grid = self.network_hwy.(linkStr).T_us_cum;
                if ~isrow(time_grid)
                    time_grid = time_grid';
                end
                if ~isrow(q_in)
                    q_in = q_in';
                end
                pbEnv.setUsFlows(time_grid,q_in);
                
                % set downstream boundary condition
                time_grid = self.network_hwy.(linkStr).T_ds_cum;
                if ~isrow(time_grid)
                    time_grid = time_grid';
                end
                if ~isrow(q_out)
                    q_out = q_out';
                end
                pbEnv.setDsFlows(time_grid,q_out);
                
                %===========================================
                % specify resolution
                nx = floor((ds_position)/dx);
                dx_res = (ds_position)/nx;
                
                x_mesh_m = 0:dx_res:ds_position;
                t_mesh_s = self.t_start:dt:self.t_end;
              
                xValues = ones(size(t_mesh_s'))*(x_mesh_m);
                tValues = t_mesh_s' * ones(size(x_mesh_m));
                
                %===========================================
                % compute Moskowitz
                result = pbEnv.explSol(tValues,xValues);
                self.M.(linkStr) = result{1}';
                
            end
            
            
        end
        
        
        %===============================================================
        function computeDensityOnGrid(self, dx)
            % This function computes the density using the M solutions
            
            % compute the density on each link
            for link = self.link_labels'
                linkStr = sprintf('link_%d', link);
                
                self.k.(linkStr) = ( self.M.(linkStr)(1:end-1,:) - ...
                                     self.M.(linkStr)(2:end,:) )/dx;                
                                 
            end
            
        end
        
        
        %===============================================================
        function M_tx = computeSolutionAtTime(self, time, link, tol)
            % This function returns the absolute vehicle ID solution at time t
            % on link.
            %   - It only computes the vehicle IDs at that time, hence no need 
            %     to compute the entire density estimation diagram, faster.
            %   - It searches shockwave intersections at that time, and the 
            %     traffic density will be aggregated and separated by the 
            %     intersectoin points. Hence less initial conditions, faster.
            %   - The shockwave intersection accuracy is up to the
            %     tolerance specified and the maximum search depth.
            %     By default, the maximum search depth is 10.
            % input: 
            %       time: float, the time that we would like to compute the
            %             solution.
            %       tol: [dt_tol; dx_tol], the accuracy tolerance in
            %            seconds and meters.
            % output:
            %       M_tx: [t, x, M]; t,x,M are all columns
            %         t, x (including the start end position)
            %         M, (the absolute vehicle ID, M(0,0) = 0)
            
            % specify the search depth and tolerance
            search_depth = 0;
            
            linkStr = sprintf('link_%d', link);

            len_link = self.network_hwy.(linkStr).para_postkm*1000;
            
            pt_found = self.searchShocksOnLine(link, [time; time], ...
                [0; len_link], [NaN; NaN], search_depth, tol);
            
            % get the time, position, and vehicle id at two bounds
            pt_us = [time, 0, self.exactSolutionAtPoints(link, time, 0)];
            pt_ds = [time, len_link, self.exactSolutionAtPoints(link, time, len_link)];
            
            % all the grid points
            M_tx = [pt_us; pt_found; pt_ds];

        end
        
        
        %===============================================================
        function plotDensityOnLink(self, link, dt, dx)
            % This function plot the solution on each link
            % simply plot in three figures to check if it is correct. 
            % The time for plotting will not be counted as the
            % computational time.
            
            % first compute the density.
            self.computeDensityOnGrid(dx);
            
            % start plotting
            linkStr = sprintf('link_%d', link);
            
            t_mesh_s = self.t_start:dt:self.t_end;
            x_mesh_m = 0:dx:self.network_hwy.(linkStr).para_postkm*1000;
            
            t_dens = t_mesh_s(1:end-1)+0.5*dt;
            x_dens = x_mesh_m(1:end-1)+0.5*dx;
            
            % ===========================================
            % transform for better color representation
            k_c_tmp = self.network_hwy.(linkStr).para_kc;
            k_m_tmp = self.network_hwy.(linkStr).para_km;
            k_trans = self.mapping(self.k.(linkStr), [0 k_c_tmp; k_c_tmp k_m_tmp],...
                [0 0.5*k_m_tmp; 0.5*k_m_tmp k_m_tmp]);
            
            clims = [0, k_m_tmp];
            scrsz = get(0,'ScreenSize');
            figure('Position',[1 1 scrsz(3) scrsz(4)]);
            
            % imagesc(flipud(k_trans), clims, 'CDataMapping','scaled')
            imagesc(t_dens, x_dens, flipud(k_trans), clims);
            h = colorbar('YTick',[0 0.5*k_m_tmp k_m_tmp],...
                'YTickLabel',{'Zero density','Critical density','Max density'});
            colormap jet
            set(gca,'fontsize',20)
            title(linkStr, 'fontsize', 30);
            xlabel({'time (s)'},'fontsize',24);
            ylabel({'space (m)'},'fontsize',24);
            
            hold on
            contour(t_mesh_s, fliplr(x_mesh_m), self.M.(linkStr), 30);
            
            % reverse the y ticks
            ax = gca;
            ax.YTick = linspace(0, self.network_hwy.(linkStr).para_postkm*1000, 11);
            ax.YTickLabel = cellstr(num2str( flipud(ax.YTick(:))));

            
            
        end
        
       
    end
    
    methods (Access = private)
       
        % ===============================================================
        function [co] = columnize(~, v)
            % columnize the array v
            
            if iscolumn(v)
                co = v;
            else
                co = v';
            end
        end
        
        
        % ===============================================================
        function [k_trans] = mapping(~, k_ori, original_interval, target_interval)
            % This function maps the values from original interval to a target interval
            % input:
            %       k_ori: the original values to be mapped
            %       original_interval: nx2 matrix. Each row is an interval to be mapped
            %           to the corresponding interval in target_interval.
            %       target_interval: nx2 matrix. Each row is the target interval.
            % output:
            %       k_trans: the mapped values which has the same dimension as k_ori
            % Example:
            % k_trans = mapping(k,[0 k_c; k_c k_m],[0 0.5*k_m; 0.5*k_m k_m]);
            % mapping (0~k_c)(k_c~k_,m) ==> (0~0.5*k_m)(0.5*k_m~k_m)
            
            k_trans = 0.0*k_ori;
            for i = 1:size(original_interval,1)
                isInDomain = (k_ori>=original_interval(i,1) & k_ori<=original_interval(i,2)+10e-6);
                k_trans(isInDomain) = target_interval(i,1) + ...
                    (k_ori(isInDomain)-original_interval(i,1))*(target_interval(i,2)-target_interval(i,1))...
                    /(original_interval(i,2)-original_interval(i,1));
            end
            
            
        end
        
        
        %===============================================================
        % Construct and solve a convex optimization program
        %===============================================================
        function setConstraintsAtMerge(self)
            % This function sets the constraints of junction solver CP for
            % a merge. 
            
            if self.num_juncs >1
                error('Current version of the code only supports one junction.')
            end
            
            
            % At each junction. NOTE: only supports one junction if using 
            % this type of junction model.
            for junc = self.junc_labels'
                juncStr = sprintf('junc_%d', junc);
                
                % preallocate memory
                num_steps = length(self.network_junc.(juncStr).T);
                % decision variables are q1, q2, e
                % constraints are compatibility condition at each point on
                % each link
                self.Aineq = zeros( 6*num_steps, 3*num_steps);
                self.bineq = zeros( 6*num_steps, 1);
                row = 0;
                
                t_grid = self.network_junc.(juncStr).T_cum;
                dt = self.network_junc.(juncStr).T';
                
                % set the decision variable index
                % decision variables: [q1(1), q1(2), ..., q2(1), q2(2), ..., e(1), e(2)]
                % where e(i) = |q2(i) - Rq1(i)|
                inlinks = self.network_junc.(juncStr).inlabel;
                self.dv_index_max = 0;
                for link  = inlinks'
                    linkStr =sprintf('link_%d', link);
                    self.dv_index.(linkStr) = [self.dv_index_max+1;...
                                               self.dv_index_max+num_steps];
                    self.dv_index_max = self.dv_index_max+num_steps;
                end
                % This is the junction e varaible
                self.dv_index.(juncStr) = [self.dv_index_max+1;...
                                           self.dv_index_max+num_steps];
                self.dv_index_max = self.dv_index_max+num_steps;
                
                % Set the sending constraints
                for link = inlinks'
                    linkStr =sprintf('link_%d', link);
                    max_sending_veh = self.maximumSendingVehiclesOnLink(link, t_grid);
                    
                    % add constraints to the Aineq*x <= bineq
                    for step = 1:num_steps
                        row = row+1;
                        self.Aineq(row, self.dv_index.(linkStr)(1) :...
                                    self.dv_index.(linkStr)(1)-1 + step) = dt(1:step);
                        self.bineq(row, 1) = max_sending_veh(step+1);
                    end
                end
                
                % Set the receiving constraints
                max_receiving_veh = self.maximumReceivingVehiclesOnLink(...
                    self.network_junc.(juncStr).outlabel, t_grid);
                % add constraints to the Aineq*x <= bineq
                for step = 1:num_steps
                    row = row+1;
                    % summation of two incoming links
                    for link = inlinks'
                        linkStr =sprintf('link_%d', link);
                        self.Aineq(row, self.dv_index.(linkStr)(1) :...
                                self.dv_index.(linkStr)(1)-1 + step) = dt(1:step);
                    end
                    self.bineq(row, 1) = max_receiving_veh(step+1);
                end
                
                % set the capacity constraints
                % q1+q2 <= q3_max; (q1<=q1_max, and q2<=q2_max are set in self.ub)
                for step = 1:num_steps
                    row = row+1;
                    % summation of two incoming links
                    for link = inlinks'
                        linkStr =sprintf('link_%d', link);
                        self.Aineq(row, self.dv_index.(linkStr)(1)-1 + step) = 1;
                    end
                    outlinkStr = sprintf('link_%d', self.network_junc.(juncStr).outlabel);
                    self.bineq(row, 1) = self.network_hwy.(outlinkStr).para_qmax;
                end
                
                % set the ratio constraints:
                % e = |q2 - Pq1| => -e <= q2-Pq1 <= e, e>=0
                for step = 1:num_steps
                    row = row+1;
                    
                    P = self.network_junc.(juncStr).ratio(2)/...
                        self.network_junc.(juncStr).ratio(1);
                    
                    inlinkStr1 = sprintf('link_%d', ...
                        self.network_junc.(juncStr).inlabel(1));
                    inlinkStr2 = sprintf('link_%d', ...
                        self.network_junc.(juncStr).inlabel(2));
                    
                    % first row, q2-Pq1-e <=0
                    self.Aineq(row, self.dv_index.(inlinkStr2)(1)-1+step) = 1;
                    self.Aineq(row, self.dv_index.(inlinkStr1)(1)-1+step) = -P;
                    self.Aineq(row, self.dv_index.(juncStr)(1)-1+step) = -1;
                    self.bineq(row, 1) = 0;
                    
                    row = row+1;
                    % first row, -q2+Pq1-e <=0
                    self.Aineq(row, self.dv_index.(inlinkStr2)(1)-1+step) = -1;
                    self.Aineq(row, self.dv_index.(inlinkStr1)(1)-1+step) = P;
                    self.Aineq(row, self.dv_index.(juncStr)(1)-1+step) = -1;
                    self.bineq(row, 1) = 0;
                    
                end
                
                
            end
            
            % set the boundary and the type
            self.lb = zeros( self.dv_index_max, 1);
            self.ub = zeros( self.dv_index_max, 1);
            for link = inlinks'
                linkStr =sprintf('link_%d', link);
                self.ub( self.dv_index.(linkStr)(1): self.dv_index.(linkStr)(2) , 1) =...
                    self.network_hwy.(linkStr).para_qmax;
            end
            % set the upper bound of e unrealistically high: 100 veh/s
            self.ub( self.dv_index.(juncStr)(1): self.dv_index.(juncStr)(2) , 1) =...
                    100;
            self.ctype(1:self.dv_index_max) = 'C';
            
        end
        
        
        %===============================================================
        function applyEntropyConAtMerge(self, junction)
            % This function applies the exact condition to construct
            % the objective function.
            % input:
            %       junction: nx1 vector, the list of junction ids that
            %       need to apply exact conditions. 
            
            if strcmp(junction, 'all')
                junction = self.junc_labels;
            end
            
            for junc = junction'
                                
                juncStr = sprintf('junc_%d',junc);
                
                % parameters for setting entropic condition at junctions
                num_steps = length(self.network_junc.(juncStr).T);
                                
                % initialize the objective function
                % initialize with 0
                self.f = zeros(self.dv_index_max,1);
                
                % if merge/diverge 
                % For each Merge or diverge junction, 
                % we add varaibel e = |q2-Rq1| for each step
                % Intuition:
                % min sum  -w(i)T(i)*(alpha*(q1(i)+q2(i)) - beta*|q2(i)-R*q1(i)|))...
                % need w(i) > w(i+1) while alpha and beta satisfying conditions
                if strcmp(self.network_junc.(juncStr).type_junc, 'merge')  
                    
                    % set entropy for this junction
                    f_junc = zeros(self.dv_index_max ,1);
                    
                    % parameter conditions
                    R_priority = self.network_junc.(juncStr).ratio(2)...
                        /self.network_junc.(juncStr).ratio(1);
                    
                    beta = 1;
                    
                    % compute the alpha and exponential factor
                    % see derive_parameters.pdf for details.
                    T_junc = self.network_junc.(juncStr).T;
                    if R_priority <= 1
                        % Given beta = 1
                        alpha = 2 + R_priority;
                        E = (alpha + 1)/(1+R_priority) + 0.01;
                    else
                        % given beta = 1
                        alpha = 1 + 2*R_priority;
                        E = (alpha + R_priority)/(1+R_priority) + 0.01;
                    end

                    % generate the weight vector
                    weight = 0.01*ones(num_steps, 1);
                    for j = 1:num_steps-1
                        weight(j+1) = E*weight(j);
                    end
                    weight(1:num_steps) = weight(num_steps:-1:1);
                    
                    % add entropic objective component
                    % sum -w(i)T(i)alpha{q1(i)+q2(i)}
                    % entropy condition:
                    for link = self.network_junc.(juncStr).inlabel'
                        linkStr = sprintf('link_%d', link);
                        f_junc(self.dv_index.(linkStr)(1):...
                            self.dv_index.(linkStr)(2),1) =...
                            weight.*T_junc;
                    end
                    
                    self.f = self.f - alpha*f_junc;
                    
                    % add the component: sum w(i) x T(i) x beta x e(i)
                    % Use e = |q2-Rq1|
                    self.f = self.f + beta*[ zeros(self.dv_index.(juncStr)(1)-1,1);...
                        weight.*T_junc;...
                        zeros( self.dv_index_max - self.dv_index.(juncStr)(2) ,1)];
                
                else
                    
                    error('Error: Current toolbox release ONLY support merge type: Junction!\n')
                    
                end
                    
                  
            
            end
            
            
        end
        
        
        %===============================================================
        function [x, fval, exitflag, output] = solveProgram(self)
            
            % Solves CP using cplex
            % min x'*H*x + f*x
            % s.t. Aineq*x <= bineq
            %      Aeq*x = beq
            %      lb <= x <= ub
            
            try
                % CPLEX
                % min f
                % s.t. Aineq*x <= bineq
                %      Aeq*x = beq
                %                 options = cplexoptimset('Diagnostics', 'on', 'TolXInteger', 10e-10);
                %                 [x, fval, exitflag, output] = ...
                %                     cplexmiqp(obj.H, obj.f, obj.Aineq, obj.bineq, obj.Aeq, obj.beq,...
                %                     [ ], [ ], [ ], obj.lb, obj.ub, obj.ctype, [ ], options);
                %                 [x, fval, exitflag, output] = ...
                %                     cplexqp(obj.H, obj.f, obj.Aineq, obj.bineq, obj.Aeq, obj.beq,...
                %                      obj.lb, obj.ub);
                [x, fval, exitflag, output] = ...
                    cplexlp(self.f, self.Aineq, self.bineq, [], [],...
                    self.lb, self.ub);
                
                fprintf ('\nSolution status = %s \n', output.cplexstatusstring);
                fprintf ('Solution value = %f \n', fval);
            catch m
                throw (m);
            end
            
            if (~isempty(x))
                %fprintf('\n Objective function:\n Min:\n');
            end
            
        end
        
        
        
        % ===============================================================
        % ===============================================================
        % Explaination of the following three blocks of code. 
        % The main idea is to utilize the lax-hopf formula and the
        % inf-morphism property to compute the exact solutions.
        % On each link, three operations may be performed:
        %   - Compute the maximum sending or receiving vehicle at the
        %   boudnary of the link, using the initial condition and the
        %   upstream OR downstream boudnary condition.
        %   - Compute the stepwise exact solution at a boundary step, 
        %   using the initial condition, the upstream and downstream
        %   boundary condition UP TO the current step.
        %   - Compute the exact solution at any point (t,x), using the
        %   initial condition, the exact upstream AND downstream boundary
        %   conditions.
        % ===============================================================
        % ===============================================================
        
        
        % ===============================================================
        % - Compute the maximum sending or receiving vehicle at the
        %   boudnary of the link, using the initial condition and the
        %   upstream OR downstream boudnary condition.
        % ===============================================================
        function rel_M = maximumSendingVehiclesOnLink(self, link, t_grid)
            % This function computes the maximum number of vehicle that can
            % be sent on link at time grid poitns
            % input: 
            %       link: the link id.
            %       t_grid: column vector, including starting and end time
            % output:
            %       rel_M: column vector same length as t_grid.
            %           cumulative relative number of vehicles that can be
            %           sent by each time grid point. rel_M(1) = 0;
            
            % save the absolute vehicle id on this link in M = [ 0; 0 ]
            rel_M = ones(length(t_grid),1)*NaN;            
            
            %===============================================================
            % extract the fundamental diagram
            linkStr = sprintf('link_%d',link);
            v_f = self.network_hwy.(linkStr).para_vf;
            w = self.network_hwy.(linkStr).para_w;    % < 0
            k_c = self.network_hwy.(linkStr).para_kc;
            k_m = self.network_hwy.(linkStr).para_km;     
            len_link = self.network_hwy.(linkStr).para_postkm*1000;
            
            % Get the initial number of vehicles 
            IC = [];
            IC(:,1) = self.network_hwy.(linkStr).X_grid_cum(1:end-1);
            IC(:,2) = self.network_hwy.(linkStr).X_grid_cum(2:end);            
            IC(:,3) = self.network_hwy.(linkStr).IC;
            
            IC_num_veh = - ( IC(:,2) - IC(:,1)).*IC(:,3);
            IC_cum_num_veh = [0; cumsum(IC_num_veh) ];
            IC_total_num_veh = IC_cum_num_veh(end);  %  should be negative
            
            % remove the last entry, which overlaps with the downstream
            % conditions
            IC_cum_num_veh(end) = [];
   
            if ~iscolumn(IC_cum_num_veh)
                IC_cum_num_veh = IC_cum_num_veh';
            end

            % sample each point
            for i = 1:length(t_grid)

                time = t_grid(i);
                position = len_link;
                
                %==========================================================
                % - Since we are computing the sending flow, hence only need
                % to consider the upstream boudnary condition
                % - It suffice to only consider the initial and upstream 
                % boundary conditions prior to the time of the being sampled point.
                num_us_pre_steps = sum(self.network_hwy.(linkStr).T_us_cum < time);
                
                if num_us_pre_steps ~= 0
                    % First extract the boundary conditions that to be considered
                    % for computing the vehilce label at the sampling point
                    BC_us = [];
                    BC_us(:,1) = self.network_hwy.(linkStr).T_us_cum(1:num_us_pre_steps);
                    BC_us(:,2) = self.network_hwy.(linkStr).T_us_cum(2:num_us_pre_steps+1);
                    % Overwrite the last time entry to the current time.
                    BC_us(num_us_pre_steps, 2) = time;
                    BC_us(:,3) = self.network_hwy.(linkStr).BC_us(1:num_us_pre_steps);

                    % Second compute the absolute vehicle ID at boundaries
                    BC_us_num_veh = (BC_us(:,2) - BC_us(:,1)).*BC_us(:,3);
                    BC_cum_us_M = [0; cumsum(BC_us_num_veh)];
                    % last entry not needed in terms of computing solutions
                    BC_cum_us_M(end) = [];
                end
                
                %==========================================================
                % compute solutions
                %==========================================================
                % compute the solution associated to the initial condition
                lowerThanKc = (IC(:,3) <= k_c);
                greaterThanKc = (IC(:,3) > k_c);
                
                %INITIAL==== <k_c ===========================================
                % first check those conditions with <= k_c
                inCharacterDomain = IC(:,1)+v_f*time <= position &...
                                    IC(:,2)+v_f*time >= position;
                
                % only applies to those <k_c
                activeCharacterDomain = lowerThanKc & inCharacterDomain;
                
                inFanDomain = IC(:,1)+v_f*time > position &...
                              IC(:,1)+w*time <= position;
                
                activeFanDomain = lowerThanKc & inFanDomain;
                
                % solution in characteristic domain
                tmp_M = min( IC_cum_num_veh(activeCharacterDomain,1) -...
                    IC(activeCharacterDomain,3).*(position - v_f*time - ...
                    IC(activeCharacterDomain,1)) );
                rel_M(i) = self.minNonEmpty(rel_M(i), tmp_M);
                % solution in fan domain
                tmp_M = min( IC_cum_num_veh(activeFanDomain,1) -...
                    k_c.*(position - v_f*time - ...
                    IC(activeFanDomain,1)  )  );
                rel_M(i) = self.minNonEmpty(rel_M(i), tmp_M);


                %INITIAL==== >k_c ===========================================
                % now check the initial condition with > k_c
                inCharacterDomain = IC(:,1)+w*time <= position &...
                                    IC(:,2)+w*time >= position; 
                activeCharacterDomain = greaterThanKc & inCharacterDomain;
                
                inFanDomain = IC(:,2)+v_f*time > position &...
                              IC(:,2)+w*time <= position; 
                activeFanDomain = greaterThanKc & inFanDomain;
                
                % solution in characteristic domain
                tmp_M = min( IC_cum_num_veh(activeCharacterDomain,1) -...
                    IC(activeCharacterDomain,3).*(position - w*time - ...
                    IC(activeCharacterDomain,1)) -k_m*time*w);
                rel_M(i) = self.minNonEmpty(rel_M(i), tmp_M);
                % solution in fan domain
                tmp_M = min( IC_cum_num_veh(activeFanDomain,1)+ IC_num_veh(activeFanDomain,1)...
                    - k_c.*(position - w*time - IC(activeFanDomain,2)) ...
                    - k_m*time*w);
                rel_M(i) = self.minNonEmpty(rel_M(i), tmp_M);

                
                if num_us_pre_steps ~= 0
                    %UPSTREAM==================================================
                    % now check the upstream conditions, only apply those
                    % boundary conditions that are not NaN
                    notNaN = ~isnan(BC_us(:,3));

                    inCharacterDomain = BC_us(:,1)+ position/v_f <= time &...
                                        BC_us(:,2)+ position/v_f >= time; 
                    activeCharacterDomain = notNaN & inCharacterDomain;

                    inFanDomain = BC_us(:,2)+ position/v_f < time;
                    activeFanDomain = notNaN & inFanDomain;

                    % solution in characteristic domain
                    tmp_M = min( BC_cum_us_M(activeCharacterDomain,1) +...
                        BC_us(activeCharacterDomain,3).*(time - position/v_f - ...
                        BC_us(activeCharacterDomain,1)) );
                    rel_M(i) = self.minNonEmpty(rel_M(i), tmp_M);
                    % solution in fan domain
                    tmp_M = min( BC_cum_us_M(activeFanDomain,1) +...
                        BC_us_num_veh(activeFanDomain,1) + ...
                        k_c*v_f.*(time - position/v_f - ...
                        BC_us(activeFanDomain,2)) );
                    rel_M(i) = self.minNonEmpty(rel_M(i), tmp_M);
                end
            end
            
            % Got the vehicle IDs from IC and upstream BC with reference
            % point to (t,x) = (0,0)
            % convert to relative M with reference to the first time point
            rel_M = rel_M - rel_M(1);
            
        end
        
        
        %===============================================================
        function rel_M = maximumReceivingVehiclesOnLink(self, link, t_grid)
            % This function computes the maximum number of vehicle that can
            % be received on link at time grid poitns
            % input: 
            %       link: the link id.
            %       t_grid: column vector, including starting and end time
            % output:
            %       rel_M: column vector same length as t_grid.
            %           cumulative relative number of vehicles that can be
            %           sent by each time grid point. rel_M(1) = 0;
            
            % save the absolute vehicle id on this link in M = [ 0; 0 ]
            rel_M = ones(length(t_grid),1)*NaN;            
            
            %===============================================================
            % extract the fundamental diagram
            linkStr = sprintf('link_%d',link);
            v_f = self.network_hwy.(linkStr).para_vf;
            w = self.network_hwy.(linkStr).para_w;    % < 0
            k_c = self.network_hwy.(linkStr).para_kc;
            k_m = self.network_hwy.(linkStr).para_km;     
            len_link = self.network_hwy.(linkStr).para_postkm*1000;
            
            % Get the initial number of vehicles 
            IC = [];
            IC(:,1) = self.network_hwy.(linkStr).X_grid_cum(1:end-1);
            IC(:,2) = self.network_hwy.(linkStr).X_grid_cum(2:end);
            IC(:,3) = self.network_hwy.(linkStr).IC;
            
            IC_num_veh = - ( IC(:,2) - IC(:,1)).*IC(:,3);
            IC_cum_num_veh = [0; cumsum(IC_num_veh) ];
            IC_total_num_veh = IC_cum_num_veh(end); %  should be negative
            
            % remove the last entry, which overlaps with the downstream
            % conditions
            IC_cum_num_veh(end) = [];
   
            if ~iscolumn(IC_cum_num_veh)
                IC_cum_num_veh = IC_cum_num_veh';
            end

            % sample each point
            for i = 1:length(t_grid)

                time = t_grid(i);
                position = 0;

                %==========================================================
                % - Since we are computing the sending flow, hence only need
                % to consider the upstream boudnary condition
                % - It suffice to only consider the initial and upstream 
                % boundary conditions prior to the time of the being sampled point.
                num_ds_pre_steps = sum(self.network_hwy.(linkStr).T_ds_cum < time);
                
                if num_ds_pre_steps ~= 0
                    % First extract the boundary conditions that to be considered
                    % for computing the vehilce label at the sampling point
                    BC_ds = [];
                    BC_ds(:,1) = self.network_hwy.(linkStr).T_ds_cum(1:num_ds_pre_steps);
                    BC_ds(:,2) = self.network_hwy.(linkStr).T_ds_cum(2:num_ds_pre_steps+1);
                    BC_ds(num_ds_pre_steps, 2) = time;
                    BC_ds(:,3) = self.network_hwy.(linkStr).BC_ds(1:num_ds_pre_steps);

                    % Second compute the absolute vehicle ID at boundaries
                    BC_ds_num_veh = (BC_ds(:,2) - BC_ds(:,1)).*BC_ds(:,3);
                    BC_cum_ds_M = [0; cumsum(BC_ds_num_veh) ] + IC_total_num_veh;
                    BC_cum_ds_M(end) = [];
                end

                %==========================================================
                % compute the solution associated to the initial condition
                lowerThanKc = (IC(:,3) <= k_c);
                greaterThanKc = (IC(:,3) > k_c);
                
                %INITIAL==== <k_c ===========================================
                % first check those conditions with <= k_c
                inCharacterDomain = IC(:,1)+v_f*time <= position &...
                                    IC(:,2)+v_f*time >= position;
                
                % only applies to those <k_c
                activeCharacterDomain = lowerThanKc & inCharacterDomain;
                
                inFanDomain = IC(:,1)+v_f*time > position &...
                              IC(:,1)+w*time <= position;
                
                activeFanDomain = lowerThanKc & inFanDomain;
                
                % solution in characteristic domain
                tmp_M = min( IC_cum_num_veh(activeCharacterDomain,1) -...
                    IC(activeCharacterDomain,3).*(position - v_f*time - ...
                    IC(activeCharacterDomain,1)) );
                rel_M(i) = self.minNonEmpty(rel_M(i), tmp_M);
                % solution in fan domain
                tmp_M = min( IC_cum_num_veh(activeFanDomain,1) -...
                    k_c.*(position - v_f*time - ...
                    IC(activeFanDomain,1)  )  );
                rel_M(i) = self.minNonEmpty(rel_M(i), tmp_M);


                %INITIAL==== >k_c ===========================================
                % now check the initial condition with > k_c
                inCharacterDomain = IC(:,1)+w*time <= position &...
                                    IC(:,2)+w*time >= position; 
                activeCharacterDomain = greaterThanKc & inCharacterDomain;
                
                inFanDomain = IC(:,2)+v_f*time > position &...
                              IC(:,2)+w*time <= position; 
                activeFanDomain = greaterThanKc & inFanDomain;
                
                % solution in characteristic domain
                tmp_M = min( IC_cum_num_veh(activeCharacterDomain,1) -...
                    IC(activeCharacterDomain,3).*(position - w*time - ...
                    IC(activeCharacterDomain,1)) -k_m*time*w);
                rel_M(i) = self.minNonEmpty(rel_M(i), tmp_M);
                % solution in fan domain
                tmp_M = min( IC_cum_num_veh(activeFanDomain,1)+ IC_num_veh(activeFanDomain,1)...
                    - k_c.*(position - w*time - IC(activeFanDomain,2)) ...
                    - k_m*time*w);
                rel_M(i) = self.minNonEmpty(rel_M(i), tmp_M);
                
                
                if num_ds_pre_steps ~= 0
                    %DOWNSTREAM==================================================
                    % now check the upstream conditions
                    notNaN = ~isnan(BC_ds(:,3));

                    inCharacterDomain = BC_ds(:,1) + (position - len_link)/w <= time &...
                                        BC_ds(:,2)+ (position - len_link)/w >= time;
                    activeCharacterDomain = notNaN & inCharacterDomain;

                    inFanDomain = BC_ds(:,2)+ (position - len_link)/w < time;
                    activeFanDomain = notNaN & inFanDomain;

                    % solution in characteristic domain
                    tmp_M = min( BC_cum_ds_M(activeCharacterDomain,1) +...
                        BC_ds(activeCharacterDomain,3).*(time -...
                        (position - len_link)/w - ...
                        BC_ds(activeCharacterDomain,1))...
                        - k_m*(position - len_link));
                    rel_M(i) = self.minNonEmpty(rel_M(i), tmp_M);
                    % solution in fan domain
                    tmp_M = min( BC_cum_ds_M(activeFanDomain,1) +...
                        BC_ds_num_veh(activeFanDomain,1) + ...
                        k_c*v_f.*(time - (position-len_link)/v_f - ...
                        BC_ds(activeFanDomain,2)) );

                    rel_M(i) = self.minNonEmpty(rel_M(i), tmp_M);
                end

            end
            
            % Got the vehicle IDs from IC and upstream BC with reference
            % point to (t,x) = (0,0)
            % convert to relative M with reference to the first time point
            rel_M = rel_M - rel_M(1);
        
        end
 
        
        % ===============================================================
        % - Compute the stepwise exact solution at a boundary step, 
        %   using the initial condition, the upstream and downstream
        %   boundary condition UP TO the current step.
        %===============================================================
        function d_M = maxSendingAtStep(self, t_sample, link)
            % This function cmoputes the amount of vehicle that can be sent 
            % on link between time self.t_ref and t_sample
            %   - The closest boundary condition (solution from x) was removed 
            %       since it may not be exact.
            %   - Do NOT call externally.
            % input: 
            %        t_sample: float, the boundary time to be sampled
            %        link: int, the link ID to be sampled
            % output: d_M: the amount of vehicles that can be sent
            
            % since t_ref is the reference point, so simply 0
            if t_sample == self.t_ref
                d_M = 0;
                return
            end
            
            % save the absolute vehicle id on this link in M = [ 0; 0 ]
            abs_M = ones(2,1)*NaN;
            
            % compute the M of t on this link
            %===============================================================
            % extract the fundamental diagram
            linkStr = sprintf('link_%d',link);
            v_f = self.network_hwy.(linkStr).para_vf;
            w = self.network_hwy.(linkStr).para_w;    % < 0
            k_c = self.network_hwy.(linkStr).para_kc;
            k_m = self.network_hwy.(linkStr).para_km;            
            
            % Get the initial number of vehicles 
            IC = [];
            IC(:,1) = self.network_hwy.(linkStr).X_grid_cum(1:end-1);
            IC(:,2) = self.network_hwy.(linkStr).X_grid_cum(2:end);
            IC(:,3) = self.network_hwy.(linkStr).IC;
            IC_num_veh = - ( IC(:,2) - IC(:,1)).*IC(:,3);
            IC_cum_num_veh = [0; cumsum(IC_num_veh) ];
            IC_total_num_veh = IC_cum_num_veh(end); %  should be negative
            
            % remove the last entry, which overlaps with the downstream
            % conditions
            IC_cum_num_veh(end) = [];
   
            if ~iscolumn(IC_cum_num_veh)
                IC_cum_num_veh = IC_cum_num_veh';
            end
            
            % select the boudnary conditions we need to use to compute the
            % vehicle id at the sampling point.
            % We only need to use the boundary conditions both upstream and
            % downstream up to the sampling point time.
            % the number of boundary steps we need to take into account
            num_us_pre_steps = sum(self.network_hwy.(linkStr).T_us_cum < t_sample);
            num_ds_pre_steps = sum(self.network_hwy.(linkStr).T_ds_cum < t_sample);
            
            
            % First extract the boundary conditions that to be considered
            % for computing the vehilce label at the sampling point
            BC_us = [];
            BC_us(:,1) = self.network_hwy.(linkStr).T_us_cum(1:num_us_pre_steps);
            BC_us(:,2) = self.network_hwy.(linkStr).T_us_cum(2:num_us_pre_steps+1);
            BC_us(num_us_pre_steps, 2) = t_sample;
            BC_us(:,3) = self.network_hwy.(linkStr).BC_us(1:num_us_pre_steps);
            
            BC_ds = [];
            BC_ds(:,1) = self.network_hwy.(linkStr).T_ds_cum(1:num_ds_pre_steps);
            BC_ds(:,2) = self.network_hwy.(linkStr).T_ds_cum(2:num_ds_pre_steps+1);
            BC_ds(num_ds_pre_steps, 2) = t_sample;
            BC_ds(:,3) = self.x( self.dv_index.(linkStr)(1):...
                                self.dv_index.(linkStr)(1) - 1 + num_ds_pre_steps);
            
            % Second compute the absolute vehicle ID
            BC_us_num_veh = (BC_us(:,2) - BC_us(:,1)).*BC_us(:,3);
            BC_cum_us_M = [0; cumsum(BC_us_num_veh)];
            BC_cum_us_M(end) = [];
            
            BC_ds_num_veh = (BC_ds(:,2) - BC_ds(:,1)).*BC_ds(:,3);
            BC_cum_ds_M = [0; cumsum(BC_ds_num_veh) ] + IC_total_num_veh;
            BC_cum_ds_M(end) = [];
            
            % remove the nonexact boundary condition
            if sum( self.network_hwy.(linkStr).T_ds_cum > self.t_ref &...
                    self.network_hwy.(linkStr).T_ds_cum < t_sample ) ~= 0
                % remove the last two pieces of downstream conditions
                BC_ds( max(1,num_ds_pre_steps-1):num_ds_pre_steps , 3) = NaN;
            else
                % remove the last piece of downstream condition
                BC_ds( num_ds_pre_steps , 3) = NaN;
            end
            
            % points coordinates
            len_link = self.network_hwy.(linkStr).para_postkm*1000;
            position = len_link;
            
            t_array = [self.t_ref; t_sample];
            for i = 1:length(t_array)
                %==========================================================
                % compute the solution associated to the initial condition
                lowerThanKc = (IC(:,3) <= k_c);
                greaterThanKc = (IC(:,3) > k_c);
                
                %INITIAL==== <k_c ===========================================
                % first check those conditions with <= k_c
                inCharacterDomain = IC(:,1)+v_f*t_array(i) <= position &...
                                    IC(:,2)+v_f*t_array(i) >= position;
                
                % only applies to those <k_c
                activeCharacterDomain = lowerThanKc & inCharacterDomain;
                
                inFanDomain = IC(:,1)+v_f*t_array(i) > position &...
                              IC(:,1)+w*t_array(i) <= position;
                
                activeFanDomain = lowerThanKc & inFanDomain;
                
                % solution in characteristic domain
                tmp_M = min( IC_cum_num_veh(activeCharacterDomain,1) -...
                    IC(activeCharacterDomain,3).*(position - v_f*t_array(i) - ...
                    IC(activeCharacterDomain,1)) );
                abs_M(i) = self.minNonEmpty(abs_M(i), tmp_M);
                % solution in fan domain
                tmp_M = min( IC_cum_num_veh(activeFanDomain,1) -...
                    k_c.*(position - v_f*t_array(i) - ...
                    IC(activeFanDomain,1)  )  );
                abs_M(i) = self.minNonEmpty(abs_M(i), tmp_M);
                
                %INITIAL==== >k_c ===========================================
                % now check the initial condition with > k_c
                inCharacterDomain = IC(:,1)+w*t_array(i) <= position &...
                                    IC(:,2)+w*t_array(i) >= position; 
                activeCharacterDomain = greaterThanKc & inCharacterDomain;
                
                inFanDomain = IC(:,2)+v_f*t_array(i) > position &...
                              IC(:,2)+w*t_array(i) <= position; 
                activeFanDomain = greaterThanKc & inFanDomain;
                
                % solution in characteristic domain
                tmp_M = min( IC_cum_num_veh(activeCharacterDomain,1) -...
                    IC(activeCharacterDomain,3).*(position - w*t_array(i) - ...
                    IC(activeCharacterDomain,1)) -k_m*t_array(i)*w);
                abs_M(i) = self.minNonEmpty(abs_M(i), tmp_M);
                % solution in fan domain
                tmp_M = min( IC_cum_num_veh(activeFanDomain,1)+ IC_num_veh(activeFanDomain,1)...
                    - k_c.*(position - w*t_array(i) - IC(activeFanDomain,2)) ...
                    - k_m*t_array(i)*w);
                abs_M(i) = self.minNonEmpty(abs_M(i), tmp_M);
                
                %UPSTREAM==================================================
                % now check the upstream conditions, only apply those
                % boundary conditions that are not NaN
                notNaN = ~isnan(BC_us(:,3));
                
                inCharacterDomain = BC_us(:,1)+ position/v_f <= t_array(i) &...
                                    BC_us(:,2)+ position/v_f >= t_array(i); 
                activeCharacterDomain = notNaN & inCharacterDomain;
                
                inFanDomain = BC_us(:,2)+ position/v_f < t_array(i);
                activeFanDomain = notNaN & inFanDomain;
                
                % solution in characteristic domain
                tmp_M = min( BC_cum_us_M(activeCharacterDomain,1) +...
                    BC_us(activeCharacterDomain,3).*(t_array(i) - position/v_f - ...
                    BC_us(activeCharacterDomain,1)) );
                abs_M(i) = self.minNonEmpty(abs_M(i), tmp_M);
                % solution in fan domain
                tmp_M = min( BC_cum_us_M(activeFanDomain,1) +...
                    BC_us_num_veh(activeFanDomain,1) + ...
                    k_c*v_f.*(t_array(i) - position/v_f - ...
                    BC_us(activeFanDomain,2)) );
                abs_M(i) = self.minNonEmpty(abs_M(i), tmp_M);
                
                
                %DOWNSTREAM==================================================
                % now check the upstream conditions
                notNaN = ~isnan(BC_ds(:,3));
                
                inCharacterDomain = BC_ds(:,1) + (position - len_link)/w <= t_array(i) &...
                                    BC_ds(:,2)+ (position - len_link)/w >= t_array(i);
                activeCharacterDomain = notNaN & inCharacterDomain;
                
                inFanDomain = BC_ds(:,2)+ (position - len_link)/w < t_array(i);
                activeFanDomain = notNaN & inFanDomain;
                
                % solution in characteristic domain
                tmp_M = min( BC_cum_ds_M(activeCharacterDomain,1) +...
                    BC_ds(activeCharacterDomain,3).*(t_array(i) -...
                    (position - len_link)/w - ...
                    BC_ds(activeCharacterDomain,1))...
                    - k_m*(position - len_link));
                abs_M(i) = self.minNonEmpty(abs_M(i), tmp_M);
                % solution in fan domain
                tmp_M = min( BC_cum_ds_M(activeFanDomain,1) +...
                    BC_ds_num_veh(activeFanDomain,1) + ...
                    k_c*v_f.*(t_array(i) - (position-len_link)/v_f - ...
                    BC_ds(activeFanDomain,2)) );
                
                abs_M(i) = self.minNonEmpty(abs_M(i), tmp_M);
                
                
            end % end for each t (here we just have (t_ref and t) )
            
            if all(~isnan(abs_M)) && all(~isempty(abs_M))
                
                % use car ID at t_ref as a reference
                d_M = abs_M(2) - abs_M(1);
                
            else
                error('failed to sample points');
            end
            
            
        end 
        
        
        %===============================================================
        function d_M = maxReceivingAtStep(self, t_sample, link)
            % This function cmoputes the amount of vehicle that can be 
            %   received on link between time self.t_ref and t_sample
            %   - The closest boundary condition (solution from x) was removed 
            %       since it may not be exact.
            %   - Do NOT call externally.
            % input: 
            %        t_sample: float, the boundary time to be sampled
            %        link: int, the link ID to be sampled
            % output: d_M, the amount of vehicles that can be received
            
            % since t_ref is the reference point, so simply 0
            if t_sample == self.t_ref
                d_M = 0;
                return
            end
            
            % save the absolute vehicle id on this link in M = [ 0; 0 ]
            abs_M = ones(2,1)*NaN;
            
            % compute the M of t on this link
            %===============================================================
            % extract the fundamental diagram
            linkStr = sprintf('link_%d',link);
            v_f = self.network_hwy.(linkStr).para_vf;
            w = self.network_hwy.(linkStr).para_w;    % < 0
            k_c = self.network_hwy.(linkStr).para_kc;
            k_m = self.network_hwy.(linkStr).para_km;            
            
            % Get the initial number of vehicles 
            IC = [];
            IC(:,1) = self.network_hwy.(linkStr).X_grid_cum(1:end-1);
            IC(:,2) = self.network_hwy.(linkStr).X_grid_cum(2:end);
            IC(:,3) = self.network_hwy.(linkStr).IC;
            IC_num_veh = - ( IC(:,2) - IC(:,1)).*IC(:,3);
            IC_cum_num_veh = [0; cumsum(IC_num_veh) ];
            IC_total_num_veh = IC_cum_num_veh(end); %  should be negative
            
            % remove the last entry, which overlaps with the downstream
            % conditions
            IC_cum_num_veh(end) = [];
   
            if ~iscolumn(IC_cum_num_veh)
                IC_cum_num_veh = IC_cum_num_veh';
            end
            
            % select the boudnary conditions we need to use to compute the
            % vehicle id at the sampling point.
            % We only need to use the boundary conditions both upstream and
            % downstream up to the sampling point time.
            % the number of boundary steps we need to take into account
            num_us_pre_steps = sum(self.network_hwy.(linkStr).T_us_cum < t_sample);
            num_ds_pre_steps = sum(self.network_hwy.(linkStr).T_ds_cum < t_sample);
            
            
            % First extract the boundary conditions that to be considered
            % for computing the vehilce label at the sampling point
            BC_us = [];
            BC_us(:,1) = self.network_hwy.(linkStr).T_us_cum(1:num_us_pre_steps);
            BC_us(:,2) = self.network_hwy.(linkStr).T_us_cum(2:num_us_pre_steps+1);
            BC_us(num_us_pre_steps, 2) = t_sample;
            
            % Make sure this is the downstream link of a merge
            if self.num_juncs ~= 1
                error('Current version of the code only supports one junction in the network.')
            end
            juncStr = sprintf('junc_%d', self.junc_labels);
            if link ~= self.network_junc.(juncStr).outlabel
                error('maxReceivingAtStep() should only be applied to the downstream link of the merge.')
            end
            if ~strcmp(self.network_junc.(juncStr).type_junc, 'merge')
                error('Current version of the code only supports merge junction.')
            end
            % The decision variable only contains the upstream outflows
            inlinkStr1 = sprintf('link_%d', self.network_junc.(juncStr).inlabel(1));
            inlinkStr2 = sprintf('link_%d', self.network_junc.(juncStr).inlabel(2));
            
            BC_us(:,3) = self.x( self.dv_index.(inlinkStr1)(1):...
                                 self.dv_index.(inlinkStr1)(1) - 1 + num_us_pre_steps)...
                       + self.x( self.dv_index.(inlinkStr2)(1):...
                                 self.dv_index.(inlinkStr2)(1) - 1 + num_us_pre_steps);
            
            BC_ds = [];                 
            BC_ds(:,1) = self.network_hwy.(linkStr).T_ds_cum(1:num_ds_pre_steps);
            BC_ds(:,2) = self.network_hwy.(linkStr).T_ds_cum(2:num_ds_pre_steps+1);
            BC_ds(num_ds_pre_steps, 2) = t_sample;
            BC_ds(:,3) = self.network_hwy.(linkStr).BC_ds(1:num_ds_pre_steps);
            
            % Second compute the absolute vehicle ID
            BC_us_num_veh = (BC_us(:,2) - BC_us(:,1)).*BC_us(:,3);
            BC_cum_us_M = [0; cumsum(BC_us_num_veh)];
            BC_cum_us_M(end) = [];
            
            BC_ds_num_veh = (BC_ds(:,2) - BC_ds(:,1)).*BC_ds(:,3);
            BC_cum_ds_M = [0; cumsum(BC_ds_num_veh) ] + IC_total_num_veh;
            BC_cum_ds_M(end) = [];
            
            % remove nonexact upstream boundary condition
            if sum( self.network_hwy.(linkStr).T_us_cum > self.t_ref &...
                    self.network_hwy.(linkStr).T_us_cum < t_sample) ~= 0
                % remove the last two pieces of upstream conditions
                BC_us( max(1, num_us_pre_steps-1) : num_us_pre_steps, 3) = NaN;
            else
                % remove the last piece of upstream condition
                BC_us( num_us_pre_steps, 3) = NaN;
            end
            
            len_link = self.network_hwy.(linkStr).para_postkm*1000;
            position = 0;
           
            t_array = [self.t_ref; t_sample];
            for i = 1:length(t_array)
                %==========================================================
                % compute the solution associated to the initial condition
                lowerThanKc = (IC(:,3) <= k_c);
                greaterThanKc = (IC(:,3) > k_c);
                
                %INITIAL==== <k_c ===========================================
                % first check those conditions with <= k_c
                inCharacterDomain = IC(:,1)+v_f*t_array(i) <= position &...
                                    IC(:,2)+v_f*t_array(i) >= position;
                
                % only applies to those <k_c
                activeCharacterDomain = lowerThanKc & inCharacterDomain;
                
                inFanDomain = IC(:,1)+v_f*t_array(i) > position &...
                              IC(:,1)+w*t_array(i) <= position;
                
                activeFanDomain = lowerThanKc & inFanDomain;
                
                % solution in characteristic domain
                tmp_M = min( IC_cum_num_veh(activeCharacterDomain,1) -...
                    IC(activeCharacterDomain,3).*(position - v_f*t_array(i) - ...
                    IC(activeCharacterDomain,1)) );
                abs_M(i) = self.minNonEmpty(abs_M(i), tmp_M);
                % solution in fan domain
                tmp_M = min( IC_cum_num_veh(activeFanDomain,1) -...
                    k_c.*(position - v_f*t_array(i) - ...
                    IC(activeFanDomain,1)  )  );
                abs_M(i) = self.minNonEmpty(abs_M(i), tmp_M);
                
                %INITIAL==== >k_c ===========================================
                % now check the initial condition with > k_c
                inCharacterDomain = IC(:,1)+w*t_array(i) <= position &...
                                    IC(:,2)+w*t_array(i) >= position; 
                activeCharacterDomain = greaterThanKc & inCharacterDomain;
                
                inFanDomain = IC(:,2)+v_f*t_array(i) > position &...
                              IC(:,2)+w*t_array(i) <= position; 
                activeFanDomain = greaterThanKc & inFanDomain;
                
                % solution in characteristic domain
                tmp_M = min( IC_cum_num_veh(activeCharacterDomain,1) -...
                    IC(activeCharacterDomain,3).*(position - w*t_array(i) - ...
                    IC(activeCharacterDomain,1)) -k_m*t_array(i)*w);
                abs_M(i) = self.minNonEmpty(abs_M(i), tmp_M);
                % solution in fan domain
                tmp_M = min( IC_cum_num_veh(activeFanDomain,1)+ IC_num_veh(activeFanDomain,1)...
                    - k_c.*(position - w*t_array(i) - IC(activeFanDomain,2)) ...
                    - k_m*t_array(i)*w);
                abs_M(i) = self.minNonEmpty(abs_M(i), tmp_M);
                
                %UPSTREAM==================================================
                % now check the upstream conditions, only apply those
                % boundary conditions that are not NaN
                notNaN = ~isnan(BC_us(:,3));
                
                inCharacterDomain = BC_us(:,1)+ position/v_f <= t_array(i) &...
                                    BC_us(:,2)+ position/v_f >= t_array(i); 
                activeCharacterDomain = notNaN & inCharacterDomain;
                
                inFanDomain = BC_us(:,2)+ position/v_f < t_array(i);
                activeFanDomain = notNaN & inFanDomain;
                
                % solution in characteristic domain
                tmp_M = min( BC_cum_us_M(activeCharacterDomain,1) +...
                    BC_us(activeCharacterDomain,3).*(t_array(i) - position/v_f - ...
                    BC_us(activeCharacterDomain,1)) );
                abs_M(i) = self.minNonEmpty(abs_M(i), tmp_M);
                % solution in fan domain
                tmp_M = min( BC_cum_us_M(activeFanDomain,1) +...
                    BC_us_num_veh(activeFanDomain,1) + ...
                    k_c*v_f.*(t_array(i) - position/v_f - ...
                    BC_us(activeFanDomain,2)) );
                abs_M(i) = self.minNonEmpty(abs_M(i), tmp_M);
                
                
                %DOWNSTREAM==================================================
                % now check the upstream conditions
                notNaN = ~isnan(BC_ds(:,3));
                
                inCharacterDomain = BC_ds(:,1) + (position - len_link)/w <= t_array(i) &...
                                    BC_ds(:,2)+ (position - len_link)/w >= t_array(i);
                activeCharacterDomain = notNaN & inCharacterDomain;
                
                inFanDomain = BC_ds(:,2)+ (position - len_link)/w < t_array(i);
                activeFanDomain = notNaN & inFanDomain;
                
                % solution in characteristic domain
                tmp_M = min( BC_cum_ds_M(activeCharacterDomain,1) +...
                    BC_ds(activeCharacterDomain,3).*(t_array(i) -...
                    (position - len_link)/w - ...
                    BC_ds(activeCharacterDomain,1))...
                    - k_m*(position - len_link));
                abs_M(i) = self.minNonEmpty(abs_M(i), tmp_M);
                % solution in fan domain
                tmp_M = min( BC_cum_ds_M(activeFanDomain,1) +...
                    BC_ds_num_veh(activeFanDomain,1) + ...
                    k_c*v_f.*(t_array(i) - (position-len_link)/v_f - ...
                    BC_ds(activeFanDomain,2)) );
                
                abs_M(i) = self.minNonEmpty(abs_M(i), tmp_M);
                
                
            end % end for each t (here we just have (t_ref and t) )
            
            if all(~isnan(abs_M)) && all(~isempty(abs_M))
                
                % use car ID at t_ref as a reference
                d_M = abs_M(2) - abs_M(1);
                
            else
                error('failed to sample points');
            end
            
            
        end 
        
        
        %===============================================================
        function TF = onStraightLine(~, x, y)
            % This function tells whether three points are on the same straight line
            % input:
            %       x: 3x1 vector, float
            %       y: 3x1 vector, float
            % output:
            %       TF: true if three points on one straight line with small
            %           tolerance; otherwise false
            
            % if a horizontal line
            if abs(x(3)-x(2)) <= 1.0e-10 && abs(x(2)-x(1)) <= 1.0e-6
                TF = true;
                return
            end
            
            % if a vertical line
            if abs(y(3)-y(2)) <= 1.0e-10 && abs(y(2)-y(1)) <= 1.0e-6
                TF = true;
                return
            end
            
            % otherwise, compute slope
            slope = (y(3)-y(1))/(x(3)-x(1));
            if abs(slope*(x(2)-x(1)) + y(1) - y(2)) <= 1.0e-6
                TF = true;
            else
                TF = false;
            end
            
        end
        
        
        %===============================================================
        function d_M = exactSolutionAtStep(self, t_sample, junc)
            % This function computes the exact solution between the 
            % self.t_ref and t_sample
            %   - it samples the sending and receiving function at all links,
            %   - then it computes the exact solution
            %   - do NOT use externally
            % input: 
            %        t_sample: float, the time point to be sampled
            %        junc: the label of the junction to be sampled
            % output: 
            %       d_M: for connection, float, the exact car ID at t_sample
            %            for merge/diverge, 2x1 float, the exact car ID at
            %            two incoming/outgoing links, with respect to self.t_ref
            
            % t_ref is the reference point where M = 0
            if t_sample == self.t_ref
                d_M = 0;
                return
            end
            
            juncStr = sprintf('junc_%d',junc);
            
            if strcmp(self.network_junc.(juncStr).type_junc,'merge')
                
                inLinks = self.network_junc.(juncStr).inlabel;
                outLink = self.network_junc.(juncStr).outlabel;
                
                % slope of the priority with y: inLinks(2); and x:
                % inlinks(1)
                sPriority = self.network_junc.(juncStr).ratio(2)/...
                    self.network_junc.(juncStr).ratio(1);
                
                % here we extract the sending and receiving function
                M_R = self.maxReceivingAtStep(t_sample, outLink);
                M_S1 = self.maxSendingAtStep(t_sample, inLinks(1));
                M_S2 = self.maxSendingAtStep(t_sample, inLinks(2));
                
                % case 1: M_S1 + M_S2 <= M_R, then through flow is
                % M_s1 + M_s2
                if M_S1 + M_S2 <= M_R
                    d_M = [M_S1; M_S2];
                
                % case 2: intersection point in side (0-M_S1, 0-M_S2) box
                elseif ( M_R/(1+sPriority) <= M_S1) &&...
                        (sPriority*M_R/(1+sPriority) <= M_S2)
                    
                    d_M = [M_R/(1+sPriority); sPriority*M_R/(1+sPriority)];                    
                
                % case 3: constrained by M_S1 flow
                elseif M_R/(1+sPriority) > M_S1
                    d_M = [M_S1, M_R-M_S1];
                
                % case 4: constrained by M_S2 flow
                elseif sPriority*M_R/(1+sPriority) > M_S2
                    d_M = [M_R-M_S2 , M_S2];
                end

            else
                error('exactSolutionAtStep currently only supports Merge junction.')
            end
            
        end
        
        
        %===============================================================
        function [TF, steps] = checkSolution(self, exactTol)
            % Check if the solution is exact.
            %   - uses exactSolutionAtStep to compute the exact
            %       solution
            %   - uses sending and receiving functions to check if the flow
            %       is maximized while minimizing the priority ratio.
            % input: 
            %       net: the network object 
            %       exact_tol: cars, the tolerance for exact, for instance, 
            %           if the error to exact solution is only 1 car, we
            %           consider this solution as an exact solution and stop
            %           updating the discretization grid.
            % output: 
            %       TF: true or false, this solution is an exact solution?
            %       steps: struct, .(juncStr), the id of the first non exact
            %           step
                        
            TF = true;
            steps = struct;
            
            % check each junction
            for junc = self.junc_labels'
                
                juncStr = sprintf('junc_%d',junc);
                
                % initialize as empty
                steps.(juncStr) = [];
                
                % the boudnary discretization grid at this junction;
                % they are also saved in the links
                T_grid = self.network_junc.(juncStr).T;
                T_cum_grid = self.network_junc.(juncStr).T_cum;
                
                if strcmp(self.network_junc.(juncStr).type_junc,'merge')
                    
                    % extract the outflow of the two incoming links
                    inlinks = self.network_junc.(juncStr).inlabel;
                    outlink = self.network_junc.(juncStr).outlabel;
                    linkStr = sprintf('link_%d',inlinks(1));
                    q_s1 = self.x(self.dv_index.(linkStr)(1):...
                                 self.dv_index.(linkStr)(2) );
                    linkStr = sprintf('link_%d',inlinks(2));
                    q_s2 = self.x(self.dv_index.(linkStr)(1):...
                                 self.dv_index.(linkStr)(2) );
                    
                    % check each step
                    for i = 1: length(T_grid)
                        
                        self.t_ref = T_cum_grid(i);   % start time of this point
                        t_check_start = T_cum_grid(i);
                        t_check_end = T_cum_grid(i+1); % end time of this point
                       
                        % a 1 x 2 vector true unique solution from the two incoming links
                        d_M = self.exactSolutionAtStep( t_check_end, junc);
                        
                        if abs( q_s1(i)*T_grid(i)-d_M(1) ) > exactTol ||...
                                abs( q_s2(i)*T_grid(i)-d_M(2) ) > exactTol
                            % The solution is not exact caused by
                            % potentially two reasons:
                            % - There are shocks in the sending or
                            %   receiving on connecting links at this step.
                            % - There are no shocks. But the convex program
                            %   produced non-exact solution due to control
                            %   action or failed entropy condition
                            d_M_SR = [self.maxSendingAtStep(t_check_end, inlinks(1));
                                      self.maxSendingAtStep(t_check_end, inlinks(2));
                                      self.maxReceivingAtStep(t_check_end, outlink)];
                            t_C = (t_check_start+t_check_end)/2;
                            d_M_C = [self.maxSendingAtStep(t_C, inlinks(1));
                                     self.maxSendingAtStep(t_C, inlinks(2));
                                     self.maxReceivingAtStep(t_C, outlink)];
                            
                            if self.onStraightLine([t_check_start, t_C, t_check_end]',[0, d_M_C(1), d_M_SR(1)]') &&...
                                  self.onStraightLine([t_check_start, t_C, t_check_end]',[0, d_M_C(2), d_M_SR(2)]') &&...
                                  self.onStraightLine([t_check_start, t_C, t_check_end]',[0, d_M_C(3), d_M_SR(3)]')
                                % - There are no shocks. But the convex program
                                %   produced non-exact solution due to control
                                %   action or failed entropy condition
                                warning('WARNING: Step %d is not exact\n', i);
                                continue
                            else
                                % - There are shocks in the sending or
                                %   receiving on connecting links at this step.
                                TF = false;
                                steps.(juncStr) = i;
                                break
                            end
                            
                        else
                            % If the first half sends more and second half
                            % sends less, then the average may be same as
                            % the exact solution. Hence double check the
                            % center point.
                            t_C = (t_check_start+t_check_end)/2;
                            d_M_C = self.exactSolutionAtStep( t_C, junc);
                            
                            if abs( q_s1(i)*T_grid(i)/2 - d_M_C(1) ) > exactTol ||...
                                    abs( q_s2(i)*T_grid(i)/2 - d_M_C(2) ) > exactTol
                                TF = false;
                                steps.(juncStr) = i;
                                break
                            else
                                % exact solutution
                                continue
                            end
                        end
                    end     % end each step
                    
                else
                    
                    error('Error: checkSolution currently only supports merge junction.')
                    
                end     % end diverge
                
                
            end % end each junction
            
            % steps = unique(steps);
            
        end 
        
        
        %===============================================================
        function T_new_grid = updateTimeGrid(self, steps)
            % Updates the the time grid by searching the intersection points
            %   - This function uses the old versions:
            %       searchShocksAtStep(), findShockSlopeAtStep(),
            %   - Find the frontwave intersection point to eliminate the
            %       discretization error.
            % input: 
            %       steps: struct, .(juncStr), the first non-exact step. 
            %           Note, the update has to be done iteratively, 
            %           earlier steps first
            % output: 
            %       T_new_cum: struct, .(juncStr). The updated time grid. 
                        
            % recursively search uses a binary cut method to locate the
            % intersection point. 
            % Set the following parameters to limit the number of
            % iterations: stop if we think the position is roughtly good.
            searchDepth = 0;
            d_t = 1.0e-3;    % 0.001 second
            
            % the new discretization grid to be returned.
            T_new_grid = struct;
            
            % check each junction
            for junc = self.junc_labels'
                
                juncStr = sprintf('junc_%d',junc);
                
                % get the link label at this junction
                inlink = self.network_junc.(juncStr).inlabel;
                outlink = self.network_junc.(juncStr).outlabel;
                
                T_tmp_cum = self.network_junc.(juncStr).T_cum;
                end_time = self.network_junc.(juncStr).T_cum(end);
                num_steps = length(self.network_junc.(juncStr).T);
                
                % find the start and end time of the next step of the
                % non-exact step, since the non-exact could be caused 
                % by the intersection from the following step.
                if ~isempty(steps.(juncStr))
                    
                    nonexact_step = steps.(juncStr);
                    
                    self.t_ref = self.network_junc.(juncStr).T_cum(nonexact_step);
                    t_left = self.network_junc.(juncStr).T_cum(nonexact_step);
                
                    if nonexact_step < num_steps
                        t_right = self.network_junc.(juncStr).T_cum(nonexact_step + 2);
                    else
                        t_right = self.network_junc.(juncStr).T_cum(nonexact_step+1);
                    end
                end
                    
                
                if strcmp( self.network_junc.(juncStr).type_junc, 'connection')
                    
                    if ~isempty(steps.(juncStr))
                        
                        % search the intersection of shockwaves on incoming links
                        % t_found = searchShocksAtStep(obj, t_interval, M_interval, slope, ...
                        %                      link, bound, searchDepth, dt_tol)
                        t_found_s = searchShocksAtStep(self,...
                            [t_left, t_right]', [0, NaN]', [NaN, NaN]',...
                                            inlink, 'downstream', searchDepth, d_t);
                        % search the intersection of shockwaves on outgoing links
                        t_found_r = searchShocksAtStep(self,...
                            [t_left, t_right]', [0, NaN]', [NaN, NaN]',...
                                            outlink, 'upstream', searchDepth, d_t);
                        
                        t_found = [t_found_s; t_found_r];
                    else
                        t_found = [];
                    end
                    
                    
                    % Aggregate original step intervals if the flows in two
                    % consecutive steps are not changing
                    % This reduces the number of discretized intervals in
                    % the boudnary without causing discretization error
                    linkStr = sprintf('link_%d', outlink);
                    q_thru = self.x(self.dv_index.(linkStr)(1,1):...
                                   self.dv_index.(linkStr)(2,1));
                    tmp_group = self.groupSameElement(q_thru, inf);
                    T_tmp_cum = [T_tmp_cum(tmp_group(:,1));...
                                 end_time];
                    T_tmp_cum = unique(T_tmp_cum);
                    
                % if the onrampjunc. Then find the boundary intersection
                % point at the upstream freeway and downstream freeway
                elseif strcmp( self.network_junc.(juncStr).type_junc, 'onrampjunc')
                    
                    % find the upstream freeway, onramp, and the downstream
                    % freeway
                    inlinks = self.network_junc.(juncStr).inlabel;
                    outlink = self.network_junc.(juncStr).outlabel;
                    linkStr1 = sprintf('link_%d',inlinks(1));
                    
                    if strcmp(self.network_hwy.(linkStr1).para_linktype, 'freeway')
                        up_link = inlinks(1);
                    else
                        up_link = inlinks(2);
                    end
                    
                    if ~isempty(steps.(juncStr))
                        
                        % search the intersection of shockwaves on incoming links
                        % t_found = searchShocksAtStep(obj, t_interval, M_interval, slope, ...
                        %                      link, bound, searchDepth, dt_tol)
                        t_found_s = searchShocksAtStep(self, [t_left, t_right]', [0, NaN]',...
                                            [NaN, NaN]', up_link, 'downstream', searchDepth, d_t);
                        % search the intersection of shockwaves on outgoing links
                        t_found_r = searchShocksAtStep(self, [t_left, t_right]', [0, NaN]',...
                                            [NaN, NaN]', outlink, 'upstream', searchDepth, d_t);
                        
                        t_found = [t_found_s; t_found_r];
                    else
                        t_found = [];
                    end
                    
                    % Aggregate original step intervals if the flows in two
                    % consecutive steps are not changing
                    linkStr = sprintf('link_%d', inlink(1));
                    q_s1 = self.x( self.dv_index.(linkStr)(1,1):...
                                  self.dv_index.(linkStr)(2,1) );
                    linkStr = sprintf('link_%d', inlink(2));
                    q_s2 = self.x( self.dv_index.(linkStr)(1,1):...
                                  self.dv_index.(linkStr)(2,1) );
                    
                    T_tmp_cum_s1 = T_tmp_cum;
                    tmp_group = self.groupSameElement(q_s1, inf);
                    T_tmp_cum_s1 = [ T_tmp_cum_s1(tmp_group(:,1));...
                                     end_time];
                    T_tmp_cum_s2 = T_tmp_cum;
                    tmp_group = self.groupSameElement(q_s2, inf);
                    T_tmp_cum_s2 = [ T_tmp_cum_s2(tmp_group(:,1));...
                                     end_time];
                    
                    T_tmp_cum = unique([T_tmp_cum_s1; T_tmp_cum_s2]);
                    
                % if the offrampjunc. Then find the boundary intersection
                % point at the upstream freeway and downstream freeway
                elseif strcmp( self.network_junc.(juncStr).type_junc, 'offrampjunc')
                    
                    % find the upstream freeway, onramp, and the downstream
                    % freeway
                    inlink = self.network_junc.(juncStr).inlabel;
                    outlinks = self.network_junc.(juncStr).outlabel;
                    linkStr1 = sprintf('link_%d',outlinks(1));
                    
                    if strcmp(self.network_hwy.(linkStr1).para_linktype, 'freeway')
                        down_link = outlinks(1);
                    else
                        down_link = outlinks(2);
                    end
                    
                    if ~isempty(steps.(juncStr))
                        
                        % search the intersection of shockwaves on incoming links
                        % t_found = searchShocksAtStep(obj, t_interval, M_interval, slope, ...
                        %                      link, bound, searchDepth, dt_tol)
                        t_found_s = searchShocksAtStep(self,...
                            [t_left, t_right]', [0, NaN]', [NaN, NaN]',...
                                            inlink, 'downstream', searchDepth, d_t);
                        % search the intersection of shockwaves on outgoing links
                        t_found_r = searchShocksAtStep(self,...
                            [t_left, t_right]', [0, NaN]', [NaN, NaN]',...
                                            down_link, 'upstream', searchDepth, d_t);
                        
                        t_found = [t_found_s; t_found_r];
                    else
                        t_found = [];
                    end
                    
                    % Aggregate original step intervals if the flows in two
                    % consecutive steps are not changing
                    linkStr = sprintf('link_%d', outlink(1));
                    q_r1 = self.x( self.dv_index.(linkStr)(1,1):...
                                  self.dv_index.(linkStr)(2,1) );
                    linkStr = sprintf('link_%d', outlink(2));
                    q_r2 = self.x( self.dv_index.(linkStr)(1,1):...
                                  self.dv_index.(linkStr)(2,1) );
                    
                    T_tmp_cum_r1 = T_tmp_cum;
                    tmp_group = self.groupSameElement(q_r1, inf);
                    T_tmp_cum_r1 = [ T_tmp_cum_r1(tmp_group(:,1));...
                                     end_time];
                    T_tmp_cum_r2 = T_tmp_cum;
                    tmp_group = self.groupSameElement(q_r2, inf);
                    T_tmp_cum_r2 = [ T_tmp_cum_r2(tmp_group(:,1));...
                                     end_time];
                    
                    T_tmp_cum = unique([T_tmp_cum_r1; T_tmp_cum_r2]);
                    
                
                % if merge junction, then find the intersection point 
                % bewteen shock waves at the boundary on three links
                elseif strcmp( self.network_junc.(juncStr).type_junc, 'merge')
                    
                    if ~isempty(steps.(juncStr))
                    
                        % corner points found from sending link 1, 2 and
                        % receiving link 3
                        t_found_s1 = searchShocksAtStep(self,...
                            [t_left, t_right]', [0, NaN]', [NaN, NaN]',...
                                            inlink(1), 'downstream', searchDepth, d_t);
                       
                        t_found_s2 = searchShocksAtStep(self,...
                            [t_left, t_right]', [0, NaN]', [NaN, NaN]',...
                                            inlink(2), 'downstream', searchDepth, d_t);
                                        
                        t_found_r = searchShocksAtStep(self,...
                            [t_left, t_right]', [0, NaN]', [NaN, NaN]',...
                                            outlink, 'upstream', searchDepth, d_t);
                        
                        t_found = [t_found_s1; t_found_s2; t_found_r];
                    else
                        t_found = [];
                    end
                    
                    % Aggregate original step intervals if the flows in two
                    % consecutive steps are not changing
                    linkStr = sprintf('link_%d', inlink(1));
                    q_s1 = self.x( self.dv_index.(linkStr)(1,1):...
                                  self.dv_index.(linkStr)(2,1) );
                    linkStr = sprintf('link_%d', inlink(2));
                    q_s2 = self.x( self.dv_index.(linkStr)(1,1):...
                                  self.dv_index.(linkStr)(2,1) );
                    
                    T_tmp_cum_s1 = T_tmp_cum;
                    tmp_group = self.groupSameElement(q_s1, inf);
                    T_tmp_cum_s1 = [ T_tmp_cum_s1(tmp_group(:,1));...
                                     end_time];
                    T_tmp_cum_s2 = T_tmp_cum;
                    tmp_group = self.groupSameElement(q_s2, inf);
                    T_tmp_cum_s2 = [ T_tmp_cum_s2(tmp_group(:,1));...
                                     end_time];
                    
                    T_tmp_cum = unique([T_tmp_cum_s1; T_tmp_cum_s2]);
                                                  
                elseif strcmp( self.network_junc.(juncStr).type_junc, 'diverge')
                     
                    if ~isempty(steps.(juncStr))
                        
                        % corner points found from sending link 1 and
                        % receiving link 2, 3
                        t_found_s = searchShocksAtStep(self,...
                            [t_left, t_right]', [0, NaN]', [NaN, NaN]',...
                                            inlink, 'downstream', searchDepth, d_t);
                        t_found_r1 = searchShocksAtStep(self,...
                            [t_left, t_right]', [0, NaN]', [NaN, NaN]',...
                                            outlink(1), 'upstream', searchDepth, d_t);
                        t_found_r2 = searchShocksAtStep(self,...
                            [t_left, t_right]', [0, NaN]', [NaN, NaN]',...
                                            outlink(1), 'upstream', searchDepth, d_t);
                        
                        t_found = [t_found_s; t_found_r1; t_found_r2];
                    
                    else 
                        t_found = [];
                    end
                  
                    % Aggregate original step intervals if the flows in two
                    % consecutive steps are not changing
                    linkStr = sprintf('link_%d', outlink(1));
                    q_r1 = self.x( self.dv_index.(linkStr)(1,1):...
                                  self.dv_index.(linkStr)(2,1) );
                    linkStr = sprintf('link_%d', outlink(2));
                    q_r2 = self.x( self.dv_index.(linkStr)(1,1):...
                                  self.dv_index.(linkStr)(2,1) );
                    
                    T_tmp_cum_r1 = T_tmp_cum;
                    tmp_group = self.groupSameElement(q_r1, inf);
                    T_tmp_cum_r1 = [ T_tmp_cum_r1(tmp_group(:,1));...
                                     end_time];
                    T_tmp_cum_r2 = T_tmp_cum;
                    tmp_group = self.groupSameElement(q_r2, inf);
                    T_tmp_cum_r2 = [ T_tmp_cum_r2(tmp_group(:,1));...
                                     end_time];
                    
                    T_tmp_cum = unique([T_tmp_cum_r1; T_tmp_cum_r2]);
                
                else
                    t_found = [];
                end
                
                % add new found intersection point in the discretization grid
                T_tmp_cum = self.unique_tol([T_tmp_cum; t_found(t_found<end_time)], 1.0e-8);
                
                % return the duration
                T_new_grid.(juncStr).T = T_tmp_cum(2:end) - T_tmp_cum(1:end-1);
                T_new_grid.(juncStr).T_cum = T_tmp_cum;
                
                % Also update the the grid at the connecting links
                for link = self.network_junc.(juncStr).inlabel'
                        
                    linkStr = sprintf('link_%d',link);
                    T_new_grid.(linkStr).T_ds = T_tmp_cum(2:end) - T_tmp_cum(1:end-1);
                    T_new_grid.(linkStr).T_ds_cum = T_tmp_cum;
                    T_new_grid.(linkStr).BC_ds = (T_tmp_cum(2:end) - T_tmp_cum(1:end-1))*NaN;
                    
                end
                
                for link = self.network_junc.(juncStr).outlabel'
                    
                    linkStr = sprintf('link_%d',link);
                    T_new_grid.(linkStr).T_us = T_tmp_cum(2:end) - T_tmp_cum(1:end-1);
                    T_new_grid.(linkStr).T_us_cum = T_tmp_cum;
                    T_new_grid.(linkStr).BC_us = (T_tmp_cum(2:end) - T_tmp_cum(1:end-1))*NaN;
                    
                end
                
                
            end
            
        end
        
        
        %===============================================================
        function setGridAtJunc(self, junc, T_junc)
            % This function updates the grid at the junction.
            % - Update the grid at the junction property
            % - Update the grid at the boundaries of links at the junction
            % - Update the boundary condition to [NaN] num_steps x 1
            % input:
            %       junc, int, the junc id
            %       T_junc, a vector of duration of each steps.
            
            % update the junc buondary grid
            juncStr = sprintf('junc_%d', junc);
            self.network_junc.(juncStr).T = T_junc;
            self.network_junc.(juncStr).T_cum = [0; cumsum(T_junc)];
            
            % update the link boundary grid
            for link = self.network_junc.(juncStr).inlabel'
                linkStr = sprintf('link_%d', link);
                num_steps = length(T_junc);
                self.network_hwy.(linkStr).T_ds = T_junc;
                self.network_hwy.(linkStr).T_ds_cum = [0; cumsum(T_junc)];
                self.network_hwy.(linkStr).BC_ds = ones(num_steps, 1)*NaN;
            end
            
            % update the link boundary grid
            for link = self.network_junc.(juncStr).outlabel'
                linkStr = sprintf('link_%d', link);
                num_steps = length(T_junc);
                self.network_hwy.(linkStr).T_us = T_junc;
                self.network_hwy.(linkStr).T_us_cum = [0; cumsum(T_junc)];
                self.network_hwy.(linkStr).BC_us = ones(num_steps, 1)*NaN;
            end
            

        end
        
        
        %===============================================================
        function t_found = searchShocksAtStep(self, t_interval, M_interval, slope, ...
                                              link, bound, searchDepth, dt_tol)
            % Find the intersection points of the sending and receiving function at the boundary
            %   - this is an older version. The more general function is searchShocksOnLine()                              
            %   - Do NOT call externally.
            % input: 
            %        t_interval, float vector, [t_start; t_end]
            %        M_interval, float vector, [M_start; M_end], NaN if not
            %           computed
            %        link, int, the link id which defines f
            %        bound, 'upstream', 'downstream', which defines f
            %        searchDepth, int, keeps track of the search depth
            %        dr_tol, float, the minimal resolution we will investigate
            % output: found time points 1 x n column
        
            t_found = [];                  
            % if search depth is greater than 3, 
            if searchDepth > self.max_search_depth
                return
            end

            % or the time interval is shorter than dt_res, return middle point  
            if t_interval(2)-t_interval(1) <= dt_tol
                % simply return the middle point as the approximated intersection point
                t_found = sum(t_interval)/2;
                return
            end
            
            % First, compute the two boundary values and the center point
            t_L = t_interval(1);
            t_R = t_interval(2);
            t_C = (t_L+t_R)/2;
            % If value not computed yet
            if isnan(M_interval(1))
                if strcmp(bound, 'downstream')
                    M_L = self.maxSendingAtStep(t_L, link);
                elseif strcmp(bound, 'upstream')
                    M_L = self.maxReceivingAtStep(t_L, link);
                end
            else
                M_L = M_interval(1);
            end
            if isnan(M_interval(2))
                if strcmp(bound, 'downstream')
                    M_R = self.maxSendingAtStep(t_R, link);
                elseif strcmp(bound, 'upstream')
                    M_R = self.maxReceivingAtStep(t_R, link);
                end
            else
                M_R = M_interval(2);
            end
            
            if strcmp(bound, 'downstream')
                M_C = self.maxSendingAtStep(t_C, link);
            elseif strcmp(bound, 'upstream')
                M_C = self.maxReceivingAtStep(t_C, link);
            end
            
            % Second, check if on a straight line, return []
            if self.onStraightLine([t_L; t_C; t_R], [M_L; M_C; M_R])
                % double check if on a straight line. Some times, the
                % segments are symmetric and those 3 points are on the same
                % line, but the sending function is not a line
                t_Q_1 = t_L + (t_R-t_L)/4;      % first quarter
                t_Q_3 = t_L + 3*(t_R-t_L)/4;    % third quarter
                
                if strcmp(bound, 'downstream')
                    M_Q_1 = self.maxSendingAtStep(t_Q_1, link);
                    M_Q_3 = self.maxSendingAtStep(t_Q_3, link);
                elseif strcmp(bound, 'upstream')
                    M_Q_1 = self.maxReceivingAtStep(t_Q_1, link);
                    M_Q_3 = self.maxReceivingAtStep(t_Q_3, link);
                end

                if self.onStraightLine([t_L; t_Q_1; t_C], [M_L; M_Q_1; M_C]) &&...
                   self.onStraightLine([t_C; t_Q_3; t_R], [M_C; M_Q_3; M_R])
                    % sampled 5 points on the same line, then very likely
                    % to be a straight line
                    return
                end
            end
            
            % Third, find the left and right slope and compute the
            % intersection 
            if isnan(slope(1))
                s_L = self.findShockSlopeAtStep([t_L; t_C], [M_L; M_C], 'left',...
                                 link, bound, dt_tol);
            else
                s_L = slope(1);
            end
            if isnan(slope(2))
                s_R = self.findShockSlopeAtStep([t_C; t_R], [M_C; M_R], 'right',...
                            	 link, bound, dt_tol);
            else
                s_R = slope(2);
            end
            
            % intersection point
            if abs(s_L - s_R) >= 1.0e-3
                % make sure two slopes are not parallel, otherwise directly
                % split them.
                t_insct = (M_R - M_L + s_L*t_L - s_R*t_R)/(s_L - s_R);
                M_insct = s_L*(t_insct - t_L) + M_L;
                
                if t_insct <= t_R && t_insct >= t_L
                    % if the intersection is in the time step, then
                    % validate, otherwise, move on to split
                    if strcmp(bound, 'downstream')
                        M_validate = self.maxSendingAtStep(t_insct, link);
                    elseif strcmp(bound, 'upstream')
                        M_validate = self.maxReceivingAtStep(t_insct, link);
                    end
                    
                    if abs(M_insct - M_validate) <= 1e-6
                        % found the intersection
                        t_found = [t_found; t_insct];
                        return
                    end
                end
            end
            
            % Forth, recursion, split the interval if there are more than 
            % one corners based on: 
            % 1. slopes are parallel; 2. intersection outside of
            % the interval; 3. cornerpoint not valid in sample
            t_found_left = self.searchShocksAtStep([t_L; t_C], [M_L; M_C],...
                [NaN; NaN], link, bound, searchDepth+1, dt_tol);
            t_found_right = self.searchShocksAtStep([t_C; t_R], [M_C; M_R],...
                [NaN; NaN], link, bound, searchDepth+1, dt_tol);
            
            t_found = [t_found; t_found_left; t_found_right];
                        
        end
        
        
        %===============================================================
        function slope = findShockSlopeAtStep(self, t_interval, M_interval, slope_side, ...
                                   link, bound, dt_tol)
            % Find the left or right slope at the boundary point of a piecewise linear line
            %   - Do NOT call externally
            % input: 
            %        t_interval, float vector, [t_start; t_end]
            %        M_interval, float vector, [M_start; M_end], NaN if not
            %           computed
            %        slope_side, 'left','right', which side slope
            %        link, int, the link id which defines f
            %        bound, 'upstream', 'downstream', which defines f
            %        dr_tol, float, the minimal resolution we will investigate
            % output: found time points 1 x n column
                               
            % if the time interval is too small, than stop computing slope
            % since this may give large numerical error
            if t_interval(2) - t_interval(1) <= dt_tol
                warning('WARNING: The slope can not be computed\n')
                slope = NaN;
                return
            end
            
            % First, sample boundary points values and the center point
            t_L = t_interval(1);
            t_R = t_interval(2);
            t_C = (t_L+t_R)/2;
            % If value not computed yet
            if isnan(M_interval(1))
                if strcmp(bound, 'downstream')
                    M_L = self.maxSendingAtStep(t_L, link);
                elseif strcmp(bound, 'upstream')
                    M_L = self.maxReceivingAtStep(t_L, link);
                end
            else
                M_L = M_interval(1);
            end
            if isnan(M_interval(2))
                if strcmp(bound, 'downstream')
                    M_R = self.maxSendingAtStep(t_R, link);
                elseif strcmp(bound, 'upstream')
                    M_R = self.maxReceivingAtStep(t_R, link);
                end
            else
                M_R = M_interval(2);
            end
            
            if strcmp(bound, 'downstream')
                M_C = self.maxSendingAtStep(t_C, link);
            elseif strcmp(bound, 'upstream')
                M_C = self.maxReceivingAtStep(t_C, link);
            end
            
            % Second, check if three points are on the same line
            if self.onStraightLine([t_L; t_C; t_R],[M_L; M_C; M_R])
                % if true, return the slope
                slope = (M_R-M_L)/(t_R-t_L);
            else
                if strcmp(slope_side, 'left')
                    % keep the left half to find the slope
                    slope = self.findShockSlopeAtStep([t_L; t_C], [M_L; M_C], slope_side,...
                                         link, bound, dt_tol);
                elseif strcmp(slope_side, 'right')
                    % keep the right half to find the slope
                    slope = self.findShockSlopeAtStep([t_C; t_R], [M_C; M_R], slope_side,...
                                         link, bound, dt_tol);
                end
            end                                 
        end
        
        
        %===============================================================
        function indicedGroupValue = groupSameElement(~, k_ori,minlength)
            % This function is used to group the same consecutive elements in an array into blocks
            % input:
            %       k_ori: a row vector that we would like to group same elements into
            %              a group.
            %       minlength: the minimal length in each group, use inf to disable
            %              this feature.
            % output:
            %       indicedGroupValue: [starting index, end index, block value]
            % example:
            % k=[1 1 1 1 1 1 1 2 2 2 NaN NaN NaN NaN 1 1 1 1 1 2]
            % indicedGroupValue = groupSameElement(k, 5)
            % OUTPUT:
            % indicedGroupValue =
            % [1 5 1;
            %  6 7 1;
            %  8 10 2;
            %  11 14 NaN;
            %  15 19 1
            %  20 20 2];
            
            k_vec = k_ori;
            
            if (~isrow(k_vec) && ~iscolumn(k_vec) )
                error('Check k dimension\n');
            elseif (iscolumn(k_vec))
                k_vec = k_vec';
            end
            
            if (~all(isnan(k_vec)==0))
                randvalue = randn*10;
                while ~all(k_vec~=randvalue)
                    randvalue = randn*10;
                end
                k_vec(isnan(k_vec)) = randvalue;
                randvalue = floor(randvalue*100000)/100000;
            end
            
            %eliminate possible numerical error
            k_tmp = round(k_vec*100000)/100000;
            originalSize = size(k_tmp,2);
            
            i = 1;
            indicedGroupValue = [0 0 0];    %initial offset
            while (1)
                
                i=i+1;
                [values, ivalues, ik_tmp] = unique(k_tmp,'stable');
                
                % set the starting indice
                indicedGroupValue(i,1) = indicedGroupValue(i-1,2)+1;
                % set the value
                indicedGroupValue(i,3) = values(1,1);
                
                if (size(values,2) > 1) %not last
                    if (ivalues(2,1) <= minlength)
                        % in case that we prefer a minlength: e.g. each group
                        % aggregate at most 5 elements.
                        indicedGroupValue(i,2) = ivalues(2,1)-1+ indicedGroupValue(i-1,2);
                        k_tmp(:,1:ivalues(2,1)-1) = [];
                    elseif (ivalues(2,1) > minlength)
                        indicedGroupValue(i,2) = minlength + indicedGroupValue(i-1,2);
                        % remove those processed elements
                        k_tmp(:,1:minlength) = [];
                    end
                    
                    % is the last element, and smaller than minlength
                elseif(originalSize-indicedGroupValue(i,1) <= minlength)
                    indicedGroupValue(i,2) = originalSize;
                    break;
                    
                else
                    indicedGroupValue(i,2) = minlength + indicedGroupValue(i-1,2);
                    % remove those processed elements
                    k_tmp(:,1:minlength) = [];
                end
                
            end
            
            % remove the initial offset (first row)
            indicedGroupValue(1,:) = [];
            
            if (~all(isnan(k_ori)==0))
                indicedGroupValue(indicedGroupValue(:,3)==randvalue,3) = NaN;
            end
            
            
        end
        
        
        %===============================================================
        function [uniq] = unique_tol(~, array, tol)
            % This function has the same functionality as unique() but considers two values within tol as identical. 
            % input:
            %       array: the array that we would like to find the unique values.
            %       tol: the tolerance between two unique values.
            % output: 
            %       uniq: the unique values in the array
            
            uniq = unique(array);
            
            len = length(uniq);
            index = (uniq(2:len)-uniq(1:len-1))<tol;
            
            % remove the value
            if sum(index)~=0
                uniq(index) = [];
            end
            
        end
        
        
        
        
        % ===============================================================
        %   - Compute the exact solution at any point (t,x), using the
        %   initial condition, the exact upstream AND downstream boundary
        %   conditions.
        % ===============================================================
        %===============================================================
        function pt_found = searchShocksOnLine(self, link, t_interval, pos_interval, ...
                                                    M_interval, searchDepth, tol)
            % This function returns the intersection of shock waves with a line defined by the two points, not including those two points.
            %   - it can be used to find the intersection point with horizontal or vertical line segment.
            %   - View the vehicle id as a function of time and find the corner point.
            %   - If it is a (or an approximated < dt_tol) vertical line, then 
            %    view vehicle id as a function of the position, then find the corner point.
            %   - Once the time or the position of the corner point found, then find 
            %    the correspoing position or time on the line.
            % input:
            %       link: the link ID
            %       t_interval: 2x1 float, the time of two points
            %       pos_interval: 2x1 float, the position of two points in meters
            %       M_interval: 2x1 float, the absolute vehicle ID, to save computation time
            %                   set NaN if not sampled yet
            %       searchDepth: int, the maximal search depth of the recursive function
            %       tol: 2x1 float [dt_tol; dx_tol] the minimal value that we would further
            %            split and locate the intersection point
            % output:
            %       pt_found: [t, x, veh_id]
            %                 t, x, and veh_id are all column vectors. Note
            %                 the two boundaries are not included.
            
            pt_found = [];

            % if search depth is 
            if searchDepth > self.max_search_depth
                return 
            end

            % First, compute the vheicle id at two bounds
            t_L = t_interval(1);
            t_R = t_interval(2);
            t_C = (t_L+t_R)/2;
            pos_L = pos_interval(1);
            pos_R = pos_interval(2);
            pos_C = (pos_L+pos_R)/2;
            % If value not computed yet
            if isnan(M_interval(1))
                M_L = self.exactSolutionAtPoints(link, t_L, pos_L);
            else
                M_L = M_interval(1);
            end
            if isnan(M_interval(2))
                M_R = self.exactSolutionAtPoints(link, t_R, pos_R);
            else
                M_R = M_interval(2);
            end
            M_C = self.exactSolutionAtPoints(link, t_C, pos_C);
            
            % Second, check if they are on a straight line
            % if t_interval > dt_tol, then find intersection by time
            if abs(t_interval(2) - t_interval(1)) > tol(1)

                if self.onStraightLine([t_L; t_C; t_R], [M_L; M_C; M_R])

                    % double check if on a straight line. Some times, the
                    % segments are symmetric and those 3 points are on the same
                    % line, but the sending function is not a line
                    t_Q_1 = t_L + (t_R-t_L)/4;      % first quarter
                    t_Q_3 = t_L + 3*(t_R-t_L)/4;    % third quarter
                    pos_Q_1 = pos_L + (pos_R-pos_L)/4;
                    pos_Q_3 = pos_L + 3*(pos_R-pos_L)/4;
                    
                    M_Q_1 = self.exactSolutionAtPoints(link, t_Q_1, pos_Q_1);
                    M_Q_3 = self.exactSolutionAtPoints(link, t_Q_3, pos_Q_3);
                    
                    if self.onStraightLine([t_L; t_Q_1; t_C], [M_L; M_Q_1; M_C]) &&...
                       self.onStraightLine([t_C; t_Q_3; t_R], [M_C; M_Q_3; M_R])
                        % sampled 5 points on the same line, then very likely
                        % to be a straight line
                        return
                    end
                end

                % Third, if not on a straing line, then find the left and right slope
                % because the slope may be computed using time or position, we do not 
                % pass the last computed value in recursion.
                s_L = self.findShockSlope(link, [t_L; t_C], [pos_L; pos_C], ...
                                         [M_L; M_C], 'left', 'time', tol);
                s_R = self.findShockSlope(link, [t_C; t_R], [pos_C; pos_R], ...
                                         [M_C; M_R], 'right', 'time', tol);

                % compute the intersection time
                if abs(s_L - s_R) >= 1.0e-3
                    % make sure two slopes are not parallel
                    t_insct = (M_R - M_L + s_L*t_L - s_R*t_R)/(s_L - s_R);
                    M_insct = s_L*(t_insct - t_L) + M_L;
                    
                    if t_insct <= t_R && t_insct >= t_L
                        % if the intersection is in the time step, then
                        % validate, otherwise, move on to split
                        pos_insct = (pos_R-pos_L)*(t_insct-t_L)/(t_R-t_L) + pos_L;
                        M_validate = self.exactSolutionAtPoints(link, t_insct, pos_insct);
                        if abs(M_insct - M_validate) <= 1e-6
                            % found the intersection
                            pt_found = [pt_found; t_insct, pos_insct, M_insct];
                            return
                        end
                    end
                end

                % Forth, recursion, split the interval if there are more than 
                % one corners based on: 
                % 1. slopes are parallel; 2. intersection outside of
                % the interval; 3. cornerpoint not valid in sample
                pt_found_left = self.searchShocksOnLine(link, [t_L; t_C], [pos_L; pos_C],...
                                    [M_L; M_C], searchDepth+1, tol);
                pt_found_right = self.searchShocksOnLine(link, [t_C; t_R], [pos_C; pos_R],...
                                    [M_C; M_R], searchDepth+1, tol);
                pt_found = [pt_found; pt_found_left; pt_found_right];


            % otherwise find by position
            elseif abs(pos_interval(2) - pos_interval(1)) > tol(2)
                
                if self.onStraightLine([pos_L; pos_C; pos_R], [M_L; M_C; M_R])

                    % double check if on a straight line. Some times, the
                    % segments are symmetric and those 3 points are on the same
                    % line, but the sending function is not a line
                    t_Q_1 = t_L + (t_R-t_L)/4;      % first quarter
                    t_Q_3 = t_L + 3*(t_R-t_L)/4;    % third quarter
                    pos_Q_1 = pos_L + (pos_R-pos_L)/4;
                    pos_Q_3 = pos_L + 3*(pos_R-pos_L)/4;
                    
                    M_Q_1 = self.exactSolutionAtPoints(link, t_Q_1, pos_Q_1);
                    M_Q_3 = self.exactSolutionAtPoints(link, t_Q_3, pos_Q_3);
                    
                    if self.onStraightLine([pos_L; pos_Q_1; pos_C], [M_L; M_Q_1; M_C]) &&...
                       self.onStraightLine([pos_C; pos_Q_3; pos_R], [M_C; M_Q_3; M_R])
                        % sampled 5 points on the same line, then very likely
                        % to be a straight line, return []
                        return
                    end
                end

                % Third, if not on a straing line, then find the left and right slope
                % because the slope may be computed using time or position, we do not 
                % pass the last computed value in recursion.
                s_L = self.findShockSlope(link, [t_L; t_C], [pos_L; pos_C], ...
                                         [M_L; M_C], 'left', 'position', tol);
                s_R = self.findShockSlope(link, [t_C; t_R], [pos_C; pos_R], ...
                                         [M_C; M_R], 'right', 'position', tol);

                % compute the intersection time
                if abs(s_L - s_R) >= 1.0e-3
                    % make sure two slopes are not parallel
                    pos_insct = (M_R - M_L + s_L*pos_L - s_R*pos_R)/(s_L - s_R);
                    M_insct = s_L*(pos_insct - pos_L) + M_L;
                    
                    if pos_insct <= pos_R && pos_insct >= pos_L
                        % if the intersection is in the time step, then
                        % validate, otherwise, move on to split
                        t_insct = (t_R-t_L)*(pos_insct-pos_L)/(pos_R-pos_L) + t_L;
                        M_validate = self.exactSolutionAtPoints(link, t_insct, pos_insct);
                        if abs(M_insct - M_validate) <= 1e-6
                            % found the intersection
                            pt_found = [pt_found; t_insct, pos_insct, M_insct];
                            return
                        end
                    end
                end

                % Forth, recursion, split the interval if there are more than 
                % one corners based on: 
                % 1. slopes are parallel; 2. intersection outside of
                % the interval; 3. cornerpoint not valid in sample
                pt_found_left = self.searchShocksOnLine(link, [t_L; t_C], [pos_L; pos_C],...
                                    [M_L; M_C], searchDepth+1, tol);
                pt_found_right = self.searchShocksOnLine(link, [t_C; t_R], [pos_C; pos_R],...
                                    [M_C; M_R], searchDepth+1, tol);
                pt_found = [pt_found; pt_found_left; pt_found_right];


            % if the area is too small to further split for more refined search 
            % then just return the center point
            else
                
                pt_found = [sum(t_interval)/2, sum(pos_interval)/2,...
                            self.exactSolutionAtPoints(link, sum(t_interval)/2, sum(pos_interval)/2)];
                return

            end

            

        end


        %===============================================================
        function slope = findShockSlope(self, link, t_interval, pos_interval, ...
                                       M_interval, slope_side, wrt, tol)
            % This funciton returns the slope of the vehicle id with respect to time or position.
            % This is different from findShockSlopeAtStep function since findShockSlopeAtStep uses 
            % maxSendingAtStep which removes the cloese boundary condition, while this function uses 
            % exactSolutionAtPoints which get the vehicle id at a certain points
            % input:
            %       link: the link id
            %       t_interval: 2x1 float vector
            %       pos_interval: 2x1 float vector define the line
            %       M_interval: 2x1 float vector, the absolute vehicle id
            %       slope_side: 'left', the first point side; 'right', the second point side
            %       wrt: 'time', 'position', with respect to
            %       tol: [dt_tol; dx_tol], the tolerance for stopping recursion.
            
            if strcmp(wrt, 'time')

                % if the time interval is too small, then stop computing slope since 
                % it may introduce numerical error
                if t_interval(2)-t_interval(1) <= tol(1)
                    warning('WARNING: The slope can not be computed\n')
                    slope = NaN;
                    return
                end

                % First, sample boundary points values and the center point
                t_L = t_interval(1);
                t_R = t_interval(2);
                t_C = (t_L+t_R)/2;
                pos_L = pos_interval(1);
                pos_R = pos_interval(2);
                pos_C = (pos_L+pos_R)/2;
                % If value not computed yet
                if isnan(M_interval(1))
                    M_L = self.exactSolutionAtPoints(link, t_L, pos_L);
                else
                    M_L = M_interval(1);
                end
                if isnan(M_interval(2))
                    M_R = self.exactSolutionAtPoints(link, t_R, pos_R);
                else
                    M_R = M_interval(2);
                end
                
                M_C = self.exactSolutionAtPoints(link, t_C, pos_C);

                if self.onStraightLine([t_L; t_C; t_R],[M_L; M_C; M_R])
                    % if true, return the slope
                    slope = (M_R-M_L)/(t_R-t_L);
                else
                    if strcmp(slope_side, 'left')
                        % keep the left half to find the slope
                        slope = self.findShockSlope(link, [t_L; t_C], [pos_L; pos_C],...
                                         [M_L; M_C], slope_side, wrt, tol);
                    elseif strcmp(slope_side, 'right')
                        % keep the right half to find the slope
                        slope = self.findShockSlope(link, [t_C; t_R], [pos_C; pos_R],...
                                         [M_C; M_R], slope_side, wrt, tol);
                    end
                end                           



            elseif strcmp(wrt, 'position')
                
                if pos_interval(2)-pos_interval(1) <= tol(2)
                    warning('WARNING: The slope can not be computed\n')
                    slope = NaN;
                    return
                end

                % First, sample boundary points values and the center point
                t_L = t_interval(1);
                t_R = t_interval(2);
                t_C = (t_L+t_R)/2;
                pos_L = pos_interval(1);
                pos_R = pos_interval(2);
                pos_C = (pos_L+pos_R)/2;
                % If value not computed yet
                if isnan(M_interval(1))
                    M_L = self.exactSolutionAtPoints(link, t_L, pos_L);
                else
                    M_L = M_interval(1);
                end
                if isnan(M_interval(2))
                    M_R = self.exactSolutionAtPoints(link, t_R, pos_R);
                else
                    M_R = M_interval(2);
                end
                
                M_C = self.exactSolutionAtPoints(link, t_C, pos_C);

                if self.onStraightLine([pos_L; pos_C; pos_R],[M_L; M_C; M_R])
                    % if true, return the slope
                    slope = (M_R-M_L)/(pos_R-pos_L);
                else
                    if strcmp(slope_side, 'left')
                        % keep the left half to find the slope
                        slope = self.findShockSlope(link, [t_L; t_C], [pos_L; pos_C],...
                                         [M_L; M_C], slope_side, wrt, tol);
                    elseif strcmp(slope_side, 'right')
                        % keep the right half to find the slope
                        slope = self.findShockSlope(link, [t_C; t_R], [pos_C; pos_R],...
                                         [M_C; M_R], slope_side, wrt, tol);
                    end
                end      



            end

        end


        %===============================================================
        function M = exactSolutionAtPoints(self, link, times, positions)
            % This function computes the Moskowitz solution at (t,x) given 
            % initial, upstream and downstream boudnary conditions.
            % input:
            %       link: the link id to be sampled
            %       time: column vector, the time of the sampled points
            %       position: column, the position of the point in meters
            % output:
            %       M: the absolute vehicle id at the position with M=0 at left bottom
            
            % it must be points
            if length(times) ~= length(positions)
                error('ERROR: the time and position must have the same length\n')
            end

            % save the absolute vehicle id on this link in M = [ 0; 0 ]
            M = ones(length(times),1)*NaN;            
            
            %===============================================================
            % extract the fundamental diagram
            linkStr = sprintf('link_%d',link);
            v_f = self.network_hwy.(linkStr).para_vf;
            w = self.network_hwy.(linkStr).para_w;    % < 0
            k_c = self.network_hwy.(linkStr).para_kc;
            k_m = self.network_hwy.(linkStr).para_km;     
            len_link = self.network_hwy.(linkStr).para_postkm*1000;
            
            % Get the initial number of vehicles 
            IC(:,1) = self.network_hwy.(linkStr).X_grid_cum(1:end-1);
            IC(:,2) = self.network_hwy.(linkStr).X_grid_cum(2:end);
            IC(:,3) = self.network_hwy.(linkStr).IC;
            IC_num_veh = - ( IC(:,2) - IC(:,1)).*IC(:,3);
            IC_cum_num_veh = [0; cumsum(IC_num_veh) ];
            IC_total_num_veh = IC_cum_num_veh(end); %  should be negative
            
            % remove the last entry, which overlaps with the downstream
            % conditions
            IC_cum_num_veh(end) = [];
   
            if ~iscolumn(IC_cum_num_veh)
                IC_cum_num_veh = IC_cum_num_veh';
            end

            % sample each point
            for i = 1:length(times)

                time = times(i);
                position = positions(i);

                % select the boudnary conditions we need to use to compute the
                % vehicle id at the sampling point.
                % We only need to use the boundary conditions both upstream and
                % downstream up to the sampling point time.
                % the number of boundary steps we need to take into account
                num_us_pre_steps = sum(self.network_hwy.(linkStr).T_us_cum < time);
                num_ds_pre_steps = sum(self.network_hwy.(linkStr).T_ds_cum < time);

                 % First extract the boundary conditions that to be considered
                % for computing the vehilce label at the sampling point
                BC_us(:,1) = self.network_hwy.(linkStr).T_us_cum(1:num_us_pre_steps);
                BC_us(:,2) = self.network_hwy.(linkStr).T_us_cum(2:num_us_pre_steps+1);
                BC_us(num_us_pre_steps, 2) = time;
                BC_us(:,3) = self.network_hwy.(linkStr).BC_us(1:num_us_pre_steps);
                
                BC_ds(:,1) = self.network_hwy.(linkStr).T_ds_cum(1:num_ds_pre_steps);
                BC_ds(:,2) = self.network_hwy.(linkStr).T_ds_cum(2:num_ds_pre_steps+1);
                BC_ds(num_ds_pre_steps, 2) = time;
                BC_ds(:,3) = self.network_hwy.(linkStr).BC_ds(1:num_ds_pre_steps);

                % Second compute the absolute vehicle ID at boundaries
                BC_us_num_veh = (BC_us(:,2) - BC_us(:,1)).*BC_us(:,3);
                BC_cum_us_M = [0; cumsum(BC_us_num_veh)];
                BC_cum_us_M(end) = [];
                
                BC_ds_num_veh = (BC_ds(:,2) - BC_ds(:,1)).*BC_ds(:,3);
                BC_cum_ds_M = [0; cumsum(BC_ds_num_veh) ] + IC_total_num_veh;
                BC_cum_ds_M(end) = [];

                %==========================================================
                % compute the solution associated to the initial condition
                lowerThanKc = (IC(:,3) <= k_c);
                greaterThanKc = (IC(:,3) > k_c);
                
                %INITIAL==== <k_c ===========================================
                % first check those conditions with <= k_c
                inCharacterDomain = IC(:,1)+v_f*time <= position &...
                                    IC(:,2)+v_f*time >= position;
                
                % only applies to those <k_c
                activeCharacterDomain = lowerThanKc & inCharacterDomain;
                
                inFanDomain = IC(:,1)+v_f*time > position &...
                              IC(:,1)+w*time <= position;
                
                activeFanDomain = lowerThanKc & inFanDomain;
                
                % solution in characteristic domain
                tmp_M = min( IC_cum_num_veh(activeCharacterDomain,1) -...
                    IC(activeCharacterDomain,3).*(position - v_f*time - ...
                    IC(activeCharacterDomain,1)) );
                M(i) = self.minNonEmpty(M(i), tmp_M);
                % solution in fan domain
                tmp_M = min( IC_cum_num_veh(activeFanDomain,1) -...
                    k_c.*(position - v_f*time - ...
                    IC(activeFanDomain,1)  )  );
                M(i) = self.minNonEmpty(M(i), tmp_M);


                %INITIAL==== >k_c ===========================================
                % now check the initial condition with > k_c
                inCharacterDomain = IC(:,1)+w*time <= position &...
                                    IC(:,2)+w*time >= position; 
                activeCharacterDomain = greaterThanKc & inCharacterDomain;
                
                inFanDomain = IC(:,2)+v_f*time > position &...
                              IC(:,2)+w*time <= position; 
                activeFanDomain = greaterThanKc & inFanDomain;
                
                % solution in characteristic domain
                tmp_M = min( IC_cum_num_veh(activeCharacterDomain,1) -...
                    IC(activeCharacterDomain,3).*(position - w*time - ...
                    IC(activeCharacterDomain,1)) -k_m*time*w);
                M(i) = self.minNonEmpty(M(i), tmp_M);
                % solution in fan domain
                tmp_M = min( IC_cum_num_veh(activeFanDomain,1)+ IC_num_veh(activeFanDomain,1)...
                    - k_c.*(position - w*time - IC(activeFanDomain,2)) ...
                    - k_m*time*w);
                M(i) = self.minNonEmpty(M(i), tmp_M);

                %UPSTREAM==================================================
                % now check the upstream conditions, only apply those
                % boundary conditions that are not NaN
                notNaN = ~isnan(BC_us(:,3));
                
                inCharacterDomain = BC_us(:,1)+ position/v_f <= time &...
                                    BC_us(:,2)+ position/v_f >= time; 
                activeCharacterDomain = notNaN & inCharacterDomain;
                
                inFanDomain = BC_us(:,2)+ position/v_f < time;
                activeFanDomain = notNaN & inFanDomain;

                % solution in characteristic domain
                tmp_M = min( BC_cum_us_M(activeCharacterDomain,1) +...
                    BC_us(activeCharacterDomain,3).*(time - position/v_f - ...
                    BC_us(activeCharacterDomain,1)) );
                M(i) = self.minNonEmpty(M(i), tmp_M);
                % solution in fan domain
                tmp_M = min( BC_cum_us_M(activeFanDomain,1) +...
                    BC_us_num_veh(activeFanDomain,1) + ...
                    k_c*v_f.*(time - position/v_f - ...
                    BC_us(activeFanDomain,2)) );
                M(i) = self.minNonEmpty(M(i), tmp_M);

                %DOWNSTREAM==================================================
                % now check the upstream conditions
                notNaN = ~isnan(BC_ds(:,3));
                
                inCharacterDomain = BC_ds(:,1) + (position - len_link)/w <= time &...
                                    BC_ds(:,2)+ (position - len_link)/w >= time;
                activeCharacterDomain = notNaN & inCharacterDomain;
                
                inFanDomain = BC_ds(:,2)+ (position - len_link)/w < time;
                activeFanDomain = notNaN & inFanDomain;
                
                % solution in characteristic domain
                tmp_M = min( BC_cum_ds_M(activeCharacterDomain,1) +...
                    BC_ds(activeCharacterDomain,3).*(time -...
                    (position - len_link)/w - ...
                    BC_ds(activeCharacterDomain,1))...
                    - k_m*(position - len_link));
                M(i) = self.minNonEmpty(M(i), tmp_M);
                % solution in fan domain
                tmp_M = min( BC_cum_ds_M(activeFanDomain,1) +...
                    BC_ds_num_veh(activeFanDomain,1) + ...
                    k_c*v_f.*(time - (position-len_link)/v_f - ...
                    BC_ds(activeFanDomain,2)) );
                
                M(i) = self.minNonEmpty(M(i), tmp_M);


            end
        
        end
    
        
        %===============================================================
        function v = minNonEmpty(~, v1, v2)
            % The function compares the min of two values which may be NaN. 
            % If one value is NaN, it will return the other value. If both NaN, return NaN.
            
            if isempty(v1) && ~isempty(v2)
                v = min(v2);
            end
            
            if isempty(v2) && ~isempty(v1)
                v = min(v1);
            end
            
            if isempty(v1) && isempty(v2)
                v = NaN;
            end
            
            if ~isempty(v1) && ~isempty(v2)
                v = min( min(v1),min(v2) );
            end
            
        end

        

        
    end
    
        
    
    
end

    
    






