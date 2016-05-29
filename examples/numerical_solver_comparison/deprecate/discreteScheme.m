

classdef discreteScheme < handle
    % This is the class for the discretization scheme
    % Yanning Li, May 22, 2016
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
        BC;             % a struct for saving the boudnary conditions
        
        
    end
    
    methods (Access = public)
        %===============================================================
        function self = discreteScheme(t_start, t_end)
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
            
            % the discretized boundary conditions on each link
            % self.BC.linkStr.upstream = [], a array, num_steps x 1
            self.BC = struct;
            
            
        end
        
        
        %===============================================================
        function addJunc(self, junc, inlabel, outlabel, type_junc, ratio, T)
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
            %       T: num_steps x 1 numerical array; time durations for each time step
            
            % make sure no duplicated label
            if any(self.junc_labels == junc)
                error('ERROR: Duplicated junc label: %d.\n', junc)
            end
            
            juncStr = sprintf('junc_%d', junc);
            self.network_junc.(juncStr).inlabel = inlabel;
            self.network_junc.(juncStr).outlabel = outlabel;
            self.network_junc.(juncStr).type_junc = type_junc;
            self.network_junc.(juncStr).ratio = ratio;
            
            % set the time duration and the time grid
            self.network_junc.(juncStr).T = T;
            self.network_junc.(juncStr).T_cum = [0; cumsum(T)];
            
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
            %               .(IC) num_segments x 1 float
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
            % estimated. It will be set as NaN.
            % Note T_in T_out must be set. One of it may be the standard
            % discretization at the boundary of the network. The other one may
            % be the dynamically set.
            
            if ~any(self.link_labels == link)
                error('ERROR: Link %d not defined in the network.\n', link);
            end
            
            % the flow data and time must have the same length
            if (~isempty(q_in) && length(q_in)~=length(T_in)) ||...
                    (~isempty(q_out) && length(q_out)~=length(T_out))
                error('Data length must be equal to its time discretization');
            end
            
            num_steps_in = length(T_in);
            num_steps_out = length(T_out);
            
            linkStr = sprintf('link_%d',link);
            
            % save the boundary discritization grid in struct
            self.network_hwy.(linkStr).T_us = T_in;
            self.network_hwy.(linkStr).T_us_cum = [0; cumsum(T_in)];
            
            self.network_hwy.(linkStr).T_ds = T_out;
            self.network_hwy.(linkStr).T_ds_cum = [0; cumsum(T_out)];
            
            % set the boundary conditions in absolute values veh/hr
            if ~isempty(q_in)
                self.network_hwy.(linkStr).BC_us = self.columnize(q_in);
            else
                self.network_hwy.(linkStr).BC_us = ones(num_steps_in,1)*NaN;
            end
            
            if ~isempty(q_out)
                self.network_hwy.(linkStr).BC_ds = self.columnize(q_out);
            else
                self.network_hwy.(linkStr).BC_ds = ones(num_steps_out,1)*NaN;
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
        function discretizeNet(self, dt, dx)
            % This function discretize the network
            % input:
            %       dt, dx: resolution, in s, m
            
            for link = self.link_labels'
                
                linkStr = sprintf('link_%d', link); 
                
                num_cells = round( self.network_hwy.(linkStr).para_postkm*1000/dx );
                num_steps = round( (self.t_end-self.t_start)/dt );
                
                % allocate the matrix for the states. 
                self.M.(linkStr) = zeros(num_cells+1, num_steps+1);
                
                % --------------------------------------------
                % discretize the initial condition
                idx = self.network_hwy.(linkStr).X_grid_cum/dx;
                idx_start = idx(1:end-1) + 1;
                idx_end = idx(2:end);
                ini_dens = zeros(num_cells, 1);
                for j = 1:length(idx_start)
                    % set the constant density to each interval
                    ini_dens( idx_start(j):idx_end(j) ) = ...
                        self.network_hwy.(linkStr).IC(j);                   
                end
                % compute the cumulative number of vehicles in each cell
                cum_num_veh = cumsum(ini_dens*dx);
                % set in the solution matrix
                self.M.(linkStr)(2:end,1) = -cum_num_veh;
                
                % --------------------------------------------
                % discretize the upstream boundary condition
                if all(isnan(self.network_hwy.(linkStr).BC_us))
                    % connected at the junction. DO NOT set the boundary
                    % condition
                    self.BC.(linkStr).upstream = [];
                    
                else
                    idx = self.network_hwy.(linkStr).T_us_cum/dt;
                    idx_start = idx(1:end-1) + 1;
                    idx_end = idx(2:end);
                    self.BC.(linkStr).upstream = zeros(num_steps, 1);
                    for j = 1:length(idx_start)
                        % set the constant density to each interval
                        self.BC.(linkStr).upstream( idx_start(j):idx_end(j) ) = ...
                            self.network_hwy.(linkStr).BC_us(j);
                    end 
                end
                
                % --------------------------------------------
                % discretize the downstream boundary condition
                if all(isnan(self.network_hwy.(linkStr).BC_ds))
                    % connected at the junction. DO NOT set the boundary
                    % condition
                    self.BC.(linkStr).downstream = [];
                    
                else
                    idx = self.network_hwy.(linkStr).T_us_cum/dt;
                    idx_start = idx(1:end-1) + 1;
                    idx_end = idx(2:end);
                    self.BC.(linkStr).downstream = zeros(num_steps, 1);
                    for j = 1:length(idx_start)
                        % set the constant density to each interval
                        self.BC.(linkStr).downstream( idx_start(j):idx_end(j) ) = ...
                            self.network_hwy.(linkStr).BC_ds(j);
                    end 
                end
                
            end
            
            
        end
        
        
        %===============================================================
        function computeSolution(self, dt, dx)
            % This function propagates the solution on each link and save
            % the solutions the corresponding matrix.
            % input: 
            %   dt: the time resolution, seconds
            %   dx: the space resolution, meters
            
            num_steps = round( (self.t_end-self.t_start)/dt );
            
            % propagate step by step
            for step = 2:num_steps+1
            
                % -----------------------------------------------
                % Compute the flows between cells and network boundary
                % Get the densities on each link
                density = struct;
                % flows.(linkStr) = [] num_cells+1 x1
                flows = struct;
                for link = self.link_labels'
                    
                    linkStr = sprintf('link_%d', link);
                    
                    num_cells = round( self.network_hwy.(linkStr).para_postkm*1000/dx );
                    
                    density.(linkStr) = (self.M.(linkStr)( 1:end-1 , step-1)...
                        - self.M.(linkStr)( 2:end, step-1 ) )/dx;
                    
                    q_sending = self.sendingFlow(density.(linkStr), linkStr);
                    q_receiving = self.receivingFlow(density.(linkStr), linkStr);
                    
                    % first compute the flows between cells by 
                    % min(sending, receiving)
                    flows.(linkStr) = zeros(num_cells+1,1);
                    flows.(linkStr)(2:end-1,1) = min( [q_sending(1:end-1), q_receiving(2:end) ], ...
                                                 [], 2);
                                            
                    % set the network boundary data, and 
                    % temporarily set the internal boundary flow as the
                    % supply-demand
                    if ~isempty( self.BC.(linkStr).upstream )
                        % posed as weak boundary condition.
                        flows.(linkStr)(1,1) = min( self.BC.(linkStr).upstream(step-1), q_receiving(1) );
                    else
                        flows.(linkStr)(1,1) = q_receiving(1);
                    end
                    
                    if ~isempty( self.BC.(linkStr).downstream )
                        % posed as weak boundary condition.
                        flows.(linkStr)(end,1) = min( self.BC.(linkStr).downstream(step-1), q_sending(end) );
                    else
                        flows.(linkStr)(end,1) = q_sending(end);
                    end
                    
                end
                
                % -----------------------------------------------
                % Compute the flows at the junction
                for junc = self.junc_labels'
                    
                    juncStr = sprintf('junc_%d', junc);
                    
                    % get the sending and receiving flows at the junction
                    sending_flows = [];
                    receiving_flows = [];
                    
                    for outlink = self.network_junc.(juncStr).outlabel'
                        linkStr = sprintf('link_%d',outlink);
                        receiving_flows = [receiving_flows;...
                                   flows.(linkStr)(1,1)];
                    end
                    
                    for inlink = self.network_junc.(juncStr).inlabel'
                        linkStr = sprintf('link_%d', inlink);
                        sending_flows = [sending_flows; ...
                                   flows.(linkStr)(end,1)];
                    end
                    
                    [q_send, q_receive] = self.computeInternalFlow( juncStr,...
                        sending_flows, receiving_flows);
                    
                    % update the sending and receiving in flows
                    i = 1;
                    for outlink = self.network_junc.(juncStr).outlabel'
                        linkStr = sprintf('link_%d',outlink);
                        flows.(linkStr)(1,1) = q_receive(i);
                        i = i+1;
                    end
                    
                    i = 1;
                    for inlink = self.network_junc.(juncStr).inlabel'
                        linkStr = sprintf('link_%d', inlink);
                        flows.(linkStr)(end,1) = q_send(i);
                        i = i+1;
                    end
                    
                end
                
                % -----------------------------------------------
                % Update the states in the current step
                for link = self.link_labels'
                    linkStr = sprintf('link_%d', link);
                    
                    % vehicle count change
                    d_M = flows.(linkStr)*dt;
                    % update vehicle id on grid
                    self.M.(linkStr)(:,step) = self.M.(linkStr)(:,step-1) + d_M;
                    
                end
            
            end
            
        end
        
        %===============================================================
        function computeDensity(self, dx)
            % This function computes the density using the M solutions
            
            % compute the density on each link
            for link = self.link_labels'
                linkStr = sprintf('link_%d', link);
                
                self.k.(linkStr) = ( self.M.(linkStr)(1:end-1,:) - ...
                                     self.M.(linkStr)(2:end,:) )/dx;                
                                 
            end
            
        end
        
        
        %===============================================================
        function plotDensityOnLink(self, link, dt, dx)
            % This function plot the solution on each link
            % simply plot in three figures to check if it is correct. 
            % The time for plotting will not be counted as the
            % computational time.
            
            % first compute the density.
            self.computeDensity(dx);
            
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
       
        %===============================================================
        function [co] = columnize(~, v)
            % columnize the array v
            
            if iscolumn(v)
                co = v;
            else
                co = v';
            end
        end
        
        
        %===============================================================
        function flow = sendingFlow(self, density, linkStr)
            % This funciton computes the sending flow from each cell using
            % the density in each cell.
            % input: 
            %   density: a column vector
            %   linkStr: the link string which is used to find the FD
            % output:
            %   flow: num_cells x 1. 
            
            freeflow_idx = density <= self.network_hwy.(linkStr).para_kc;
            
            flow = 0*density;
            flow(freeflow_idx) = density(freeflow_idx)*...
                self.network_hwy.(linkStr).para_vf;
            flow(~freeflow_idx) = self.network_hwy.(linkStr).para_qmax;
            
        end
        
        
        %===============================================================
        function flow = receivingFlow(self, density, linkStr)
            % This funciton computes the receiving flow from each cell using
            % the density in each cell.
            % input: 
            %   density: a column vector
            %   link_str: the link string which is used to find the FD
            % output:
            %   flow: num_cells x 1
            freeflow_idx = density <= self.network_hwy.(linkStr).para_kc;
            
            flow = 0*density;
            flow(freeflow_idx) = self.network_hwy.(linkStr).para_qmax;
            flow(~freeflow_idx) = - self.network_hwy.(linkStr).para_w*...
                ( self.network_hwy.(linkStr).para_km - density(~freeflow_idx) );
            
        end
        
        
        %===============================================================
        function [q_send, q_receive] = ...
            computeInternalFlow(self, juncStr, sending_flows, receiving_flows)
            % This function computes the internal boudnary flows using the
            % junction model.
            % input:
            %   juncStr: the junction string
            %   sending_flows: the sending flows on links in the same order
            %   as .(juncStr).inlabel
            %   receiving_flows: the receiving flows on links in the same
            %   order as .(juncStr).outlabel
            
            if strcmp(self.network_junc.(juncStr).type_junc, 'merge')
                
                q_send = [0;0];
                q_receive = 0;
                
                % Scenario 1
                if sum(sending_flows) <= receiving_flows
                    
                    q_send = sending_flows;
                    q_receive = sum(q_send);
                
                else
                    % compute the intersection point
                    % q2 = q1*P
                    P = self.network_junc.(juncStr).ratio(2)/...
                        self.network_junc.(juncStr).ratio(1);
                    
                    p_intersect = [ receiving_flows/(1+P) ,...
                                    P*receiving_flows/(1+P)];                          
                    
                    % Scenario 2
                    if p_intersect(1) <= sending_flows(1) && ...
                       p_intersect(2) <= sending_flows(2)
                        q_send = p_intersect;
                        q_receive = receiving_flows;
                        
                    end
                    
                    % Scenario 3
                    if p_intersect(1) >= sending_flows(1)
                        
                        q_send(1) = sending_flows(1);
                        q_send(2) = receiving_flows - sending_flows(1);
                        q_receive = receiving_flows;
                    
                    elseif p_intersect(2) >= sending_flows(2)
                        
                        q_send(2) = sending_flows(2);
                        q_send(1) = receiving_flows - sending_flows(2);
                        q_receive = receiving_flows;
                        
                    end
                        
                
                end
                
            elseif strcmp(self.network_junc.(juncStr).type_junc, 'diverge')
                
                q_send = 0;
                q_receive = [0;0];
                
                % Scenario 1
                if sum(receiving_flows) <= sending_flows
                    
                    q_receive = receiving_flows;
                    q_send = sum(receiving_flows);
                
                else
                    % compute the intersection point
                    % q2 = q1*D
                    D = self.network_junc.(juncStr).ratio(2)/...
                        self.network_junc.(juncStr).ratio(1);
                    
                    p_intersect = [ sending_flows/(1+P) ,...
                                    P*sending_flows/(1+P)];                          
                    
                    % Scenario 2
                    if p_intersect(1) <= receiving_flows(1) && ...
                       p_intersect(2) <= receiving_flows(2)
                   
                        q_receive = p_intersect;
                        q_send = sending_flows;
                        
                    end
                    
                    % Scenario 3
                    if p_intersect(1) >= receiving_flows(1)
                        
                        q_receive(1) = receiving_flows(1);
                        q_receive(2) = sending_flows - receiving_flows(1);
                        q_send = sending_flows;
                    
                    elseif p_intersect(2) >= receiving_flows(2)
                        
                        q_receive(2) = receiving_flows(2);
                        q_receive(1) = sending_flows - receiving_flows(2);
                        q_send = sending_flows;
                        
                    end
                end
                
                
            elseif strcmp(self.network_junc.(juncStr).type_junc, 'connection')
                
                q_send = min(sending_flows, receiving_flows);
                q_receive = min(sending_flows, receiving_flows);
                
            else
                
                error('Unsupported junction type.')
                
            end
            
            
        end
        
        
        
        %===============================================================
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
        
        
    end
    
end

    
    






