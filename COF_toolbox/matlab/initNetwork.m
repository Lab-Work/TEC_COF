

classdef initNetwork < handle
    % initNetwork constructs a network object which contains the network structure, the initial and boundary conditions
    % Yanning Li, Jan 10, 2016
    % network:
    % ---- hwy: a struct, containing information for each hwy link
    % ---- junc: a struct, containing the informatino for each junc, including
    %            the time grid at the junction
    % ---- initial and boundary conditions
    % ---- other auxiliary properties
    
    properties
        
        network_junc;   % struct, .(juncStr) contains junction info
        network_hwy;    % struct, .(linkStr), link info: .IC, .X_grid, .X_grid_cum, .para_~
                        
        junc_labels;    % array, list of junction integer labels
        link_labels;    % array, list of link integer labels
        num_juncs;      % number of junctions in the network
        num_links;      % number of links in the network
        
    end
    
    methods
        %===============================================================
        function self = initNetwork(~)
            % Construct an empty net object
                        
            self.network_junc = struct;
            self.network_hwy = struct;
            
            self.junc_labels = [];
            self.link_labels = [];
            self.num_links = 0;
            self.num_juncs = 0;
            
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
        function setInternalCon(self, link, t_interval, x_interval, r_meas)
            % Set internal condition of each link
            % input:
            %       t_interval: nx2 matrix; each row [t_min_traj, t_max_traj]
            %       x_interval: nx2 matrix; each row [x_min_traj, x_max_traj]
            %       r_meas: the passing rate of each condition
            
            if ~any(self.link_labels == link)
                error('ERROR: Link %d not defined in the network\n', link);
            end
            
            if size(t_interval,1) ~= length(r_meas) || ...
                    size(x_interval,1) ~= length(r_meas)
                error('ERROR: the internal condition is not correctly defined\n')
            end
            
            linkStr = sprintf('link_%d',link);
            self.network_hwy.(linkStr).x_min_traj = x_interval(:,1);
            self.network_hwy.(linkStr).x_max_traj = x_interval(:,2);
            self.network_hwy.(linkStr).t_min_traj = t_interval(:,1);
            self.network_hwy.(linkStr).t_max_traj = t_interval(:,2);
            self.network_hwy.(linkStr).r_meas_traj = r_meas;
            
            self.network_hwy.(linkStr).v_meas_traj = ...
                (x_interval(:,2)-x_interval(:,1))./...
                (t_interval(:,2)-t_interval(:,1));

        end
        
        
        %===============================================================
        function setDensityCon(self, link, t_dens, x_interval, dens_meas)
            % Set density condition of each link
            % input:
            %       t_dens: nx1 matrix; each row [t_dens]
            %       x_interval: nx2 matrix; each row [x_min_traj, x_max_traj]
            %       dens_meas: nx1; the absolute density
            
            if ~any(self.link_labels == link)
                error('ERROR: Link %d not defined in the network\n', link);
            end
            
            if size(t_dens,1) ~= length(dens_meas) || ...
                    size(x_interval,1) ~= length(dens_meas)
                error('ERROR: the density condition is not correctly defined\n')
            end
            
            linkStr = sprintf('link_%d',link);
            self.network_hwy.(linkStr).x_min_dens = x_interval(:,1);
            self.network_hwy.(linkStr).x_max_dens = x_interval(:,2);
            self.network_hwy.(linkStr).t_dens = t_dens;
            self.network_hwy.(linkStr).dens_meas = dens_meas;

        end
        
        
        %===============================================================
        function [co] = columnize(~, v)
            % columnize the array v
            
            if iscolumn(v)
                co = v;
            else
                co = v';
            end
        end
        
     
    end
    
end

    
    






