% Yanning Li, Sep 11, 2015

% This class sets the inequality constraints for the optimization program.
% Cumunlativley, it has the following features
% 1. It is based on the setIneqConstraints_grid.m which only adds the 
%   minimum set of constarints based on heuristic knowledge. See the code 
%   description in setIneqConstraints_grid.m for details.
% 2. The variable names are changed which makes the code much more
%   readable. The original code varialble names are defined in agreement
%   with the equations in paper. 
% 3. This version supports different levels of error for the value
%   conditions.
% 4. This version supports uneven discretization of the time and space.
% 5. This version has full support for the initial condition, the boundary
%   condition, the internal condition, the density condition. Those
%   conditions are interacting. Requries quite some effort to fix.

% The index of the variables is
% [q_us, q_ds, rho_ini, L_traj, r_traj, 
%  L_dens, rho_dens, L_queue, s_queue...
%  bool_traj_min, bool_traj_max, 
%  bool_dens_min, bool_dens_max, 
%  other auxilary variables]

% The first columns are A, the last column of the inequlity matrix is b
% the signs follows Ax >= b

% Description of value conditions:
% Boundary condition: q_us, q_ds. The boundary flows at upstream and
%       downstream position
% Initial condition: rho_ini. The initial density on the road.
% Internal trajectory condition: L_traj, r_traj. Internal condition is to 
%       incorporate the trajectory data collected by probe vehicles. 
%       Each internal condition is defined by a line segment. The position 
%       of the line is defined by two points. The decision variables 
%       associated with each internal condition are the vehicle id of that 
%       line trajectory L_traj, and the passing rate r_traj. The passing 
%       rate r_traj allows vehicles to pass this traj which is realistic in 
%       multiple-lane traffic.
% Density condition: L_dens, rho_dens. The density is similar to the
%       initial condition except that it can be defined at any time. It is
%       more useful for defining the objective. For instance, minimizing the 
%       L1 norm of adjacent density conditions, which can smooth the final
%       estimates of the traffic density. 
%       Each density condition is defined by a vetical line segment. The
%       decision variables associated with each density condition are the
%       L_dens, which is the vehicle id at the x_min point; the rho_dens
%       which is the density in that line segment. 
% Soft queue limit: L, s. The label of the vehicle at (x,t) and the slack variable.
%       constriants: L <= M^{us, ds, init}; L = M^us + s; s>=0; min s^2
%       


classdef setIneqConstraints
    
    properties
        
        dv_link;     % the index of the link
                     % .upstream; .downstream; .initial; 
                     % .internal_L; .internal_r; 
                     % .density_L; .density_rho;
                     % .queue_L; .queue_s
                     % .num_step_past_us, .num_step_past_ds
        dv_max;      % keep track of the maximal index including the bools
        dv_var_max; % the maximum index for link float variables (not the bools)
        
        num_us_con;  % num of upstream boundary conditions
        num_ds_con;  % num of downstream boundary conditions
        num_initial_con;  % num of initial conditions
        
        % num of internal conditions; internal condition is defined by two
        % points in the time space domain with a velocity and rate of
        % passing, collected by probe vehicles with the same id. 
        % Let us call this more explicitely in this code as trajectory 
        % conditions from now on. For each internal condition, we need two
        % variables, the vehicle id L and the passing rate r.
        num_internal_con;  
        
        % Density conditions is defined by two points on a vertical line,
        % with the density value in between those two points.
        num_density_con;  % num of density conditions
        
        % soft queue limit points.
        % This is not a value condtion, but for posing a soft queue limit
        % Each pt associates with M and a slack variable s
        % M <= M_{us, ds, init}, and M = M_us + s, s>=0, penalize s^2
        num_soft_queue_pt;
        
        v;  % freeflow speed
        w;  % congestion slope
        k_c;    % critical density
        k_m;    % maximal density
        q_max;
        
        start_time;
        now_time;
        end_time;   % end of simulation time window
      
        %=========================================================
        % the following are for boundary conditions
        us_pos_m;     % upstream position in meters
        ds_pos_m;    % post meter of the downstream position
        
        % up and downstream time grid
        T_us;       % upstream_steps x 1
        T_us_cum;   % cumulative; upstream_steps+1 x 1
        T_ds;
        T_ds_cum;
        
        % for boundary conditions
        qin_meas;
        qout_meas;
        
        % Use absolute error expressed in percent of normal range except
        % the passing rate r.
        % q_meas-e_his*q^max < q < q_meas+e_his*q^max
        e_max;  % error 
        
        e_his;  % error for historical data; normal q_max
        e_est;  % error for the estimated density condition; normal k_c
        e_meas_flow; % error of the measurement of boundary flows; normal q_max
        e_meas_dens;    % error of the density condition value; normalk_c
        
        % For passing rate r, use relative error
        % r_meas(1-err) <= r <= r_meas(1+err)
        e_meas_traj_r;  % error of the trajectory passing rate
        
        %=========================================================
        % The following are for initial conditions
        X_grid;       % the length of each initial condition segments
        X_grid_cum;   % the grid of each initial condition separation point
        rho_ini;    % b_max x 1
        indiced_rho_ini;    
        
        
        %=========================================================
        % for internal (trajectory) conditions
        x_min_traj;
        x_max_traj;
        t_min_traj;
        t_max_traj;
        v_meas_traj;
        r_meas_traj;
        
        % for density conditions
        x_min_dens;
        x_max_dens;
        t_dens;
        dens_meas;
        
        % for soft queue points positions
        soft_queue_x;
        soft_queue_t;
        
        % the boolean variables for each point of the internal or the
        % density condition. 
        % It is an array with the length of the number of points. Each
        % entry is the number of boolean variables that need to be added. 
        % Note the number of loolean variables are computed as the log2 of the
        % number of solutions. 
        % E.g. for point (t,x) with value condition L, it associated with 4 
        %   solutions, then the number of boolean variables is 2 which gives 
        %   four combinations. This saves memory.
        nb_min_traj;
        nb_max_traj;
        nb_min_dens;
        nb_max_dens;

        % auxiliary variables
        num_aux;
        size_row;
        
        % total number of boolean variables
        num_bool;
        
    end
    
    methods
        
        % This function set the inequality constraints of one link
        % input: 
        %       para: v, w, k_c, k_m: road parameters; postm: length in meters
        %       start_time, end_time: start and end time of the simulation
        %       now_time: the current time (the timestamp of latest data
        %           mesurement). Data before this will be measurement data;
        %           after this will be historical data.
        %       Boundary_con: a struct, with fields
        %           .BC_us: column vector, float, upstream measuremnt data
        %           .BC_ds: column vector, float, downstream measuremnt data
        %           .T_us: column vector, float, the upstream time grid, same
        %               length as qin_meas
        %           .T_ds: column vector, float, downstream time grid
        %       Initial_con: a struct with fields
        %           .X_grid_cum: column vector, float, cumulative space grid
        %           .IC: column vector, float, initial density 
        %       Traj_con: a struct with fields for trajectory condition
        %       Dens_con: a struct with fields for density condition
        %       soft_queue_limit: float, the soft limit of the length of
        %           the queue in meters. Auto set every 30 s; 
        %           Set value to NaN to disable
        %       error: struct, must have a field e_default
        %              other fields corresponding to each property
        function self = setIneqConstraints(...
                para,...
                start_time, now_time, end_time,...
                Boundary_con,...
                Initial_con, Traj_con, Dens_con, ...
                soft_queue_limit,...
                errors)            
            
            % set road parameters
            self.v = para.vf;
            self.w = para.w;
            self.k_c = para.kc;
            self.k_m = para.km;
            self.q_max = self.v*self.k_c;
            
            self.start_time = start_time;
            self.now_time = now_time;
            self.end_time = end_time;
            
            % set the errors
            if isfield(errors, 'e_his')
                self.e_his = errors.e_his;
            else
                self.e_his = errors.e_default;
            end
            if isfield(errors, 'e_est')
                self.e_est = errors.e_est;
            else
                self.e_est = errors.e_default;
            end
            if isfield(errors, 'e_meas_flow')
                self.e_meas_flow = errors.e_meas_flow;
            else
                self.e_meas_flow = errors.e_default;
            end
            if isfield(errors, 'e_meas_traj_r')
                self.e_meas_traj_r = errors.e_meas_traj_r;
            else
                self.e_meas_traj_r = errors.e_default;
            end
            if isfield(errors, 'e_meas_dens')
                self.e_meas_dens = errors.e_meas_dens;
            else
                self.e_meas_dens = errors.e_default;
            end
                
            self.us_pos_m = 0;  % relative postm for each link
            self.ds_pos_m = para.postm;
            
            
            
            self.dv_max = 0;
            self.num_us_con = size(Boundary_con.BC_us,1);   % boundary conditions
            self.dv_link.upstream = [self.dv_max+1; ...
                                     self.dv_max+self.num_us_con];
            self.dv_max = self.dv_max + self.num_us_con;
            
            self.num_ds_con = size(Boundary_con.BC_ds,1);
            self.dv_link.downstream = [self.dv_max+1;...
                                       self.dv_max+self.num_ds_con];
            self.dv_max = self.dv_max + self.num_ds_con;
                                   
            self.num_initial_con = size(Initial_con.IC,1);   % initial conditions
            self.dv_link.initial = [self.dv_max+1;...
                                    self.dv_max+self.num_initial_con];
            self.dv_max = self.dv_max + self.num_initial_con;
            
            if isfield(Traj_con, 'x_min_traj')
                self.num_internal_con = length(Traj_con.x_min_traj);
                % internal L
                self.dv_link.internal_L = [self.dv_max+1;
                                       self.dv_max+self.num_internal_con];
                self.dv_max = self.dv_max + self.num_internal_con;
                % internal r
                self.dv_link.internal_r = [self.dv_max+1;
                                       self.dv_max+self.num_internal_con];
                self.dv_max = self.dv_max + self.num_internal_con;
            else
                self.num_internal_con = 0;     
                % set as none do not set index
            end
                        
            if isfield(Dens_con, 'x_min_dens')
                self.num_density_con = length(Dens_con.x_min_dens);
                % density L
                self.dv_link.density_L = [self.dv_max+1;
                                          self.dv_max+self.num_density_con];
                self.dv_max = self.dv_max + self.num_density_con;
                % density rho
                self.dv_link.density_rho = [self.dv_max+1;
                                            self.dv_max+self.num_density_con];
                self.dv_max = self.dv_max + self.num_density_con;
            else
                self.num_density_con = 0;         
            end
            
            if ~isnan(soft_queue_limit)
                
                if self.num_internal_con ~= 0 || self.num_density_con ~= 0
                    error('Current version does not support soft queue limit with internal or density condition.')
                end
                
                tmp_grid = self.now_time:30:self.end_time;
                if tmp_grid(end)~= self.end_time
                    tmp_grid = [tmp_grid'; self.end_time];
                end
                self.soft_queue_t = tmp_grid';
                self.num_soft_queue_pt = length(tmp_grid);
                self.soft_queue_x = ones(self.num_soft_queue_pt,1)*...
                    (self.ds_pos_m-soft_queue_limit);
                % soft queue L
                self.dv_link.queue_L = [self.dv_max+1;
                                        self.dv_max+self.num_soft_queue_pt];
                self.dv_max = self.dv_max + self.num_soft_queue_pt;   
                % soft queue s
                self.dv_link.queue_s = [self.dv_max+1; 
                                        self.dv_max+self.num_soft_queue_pt];
                self.dv_max = self.dv_max + self.num_soft_queue_pt;
            else
                self.num_soft_queue_pt = 0;
            end            
            self.dv_var_max = self.dv_max;
            
            
            % Boundary condition
            self.T_us = Boundary_con.T_us;
            self.T_us_cum = [0; cumsum(Boundary_con.T_us)];
            self.T_ds = Boundary_con.T_ds;
            self.T_ds_cum = [0; cumsum(Boundary_con.T_ds)];
            self.qin_meas = Boundary_con.BC_us;
            self.qout_meas = Boundary_con.BC_ds;
            
            % find out the number of past steps and the number of predict
            % steps
            self.dv_link.num_step_past_us = sum(self.T_us_cum <= self.now_time)-1;
            self.dv_link.num_step_past_ds = sum(self.T_ds_cum <= self.now_time)-1;
            
            % internal trajectory condition
            if isempty(Traj_con) || isempty(fieldnames(Traj_con))
                % if the trajectory condition struct is empty, then disable
                % trajectory condition
                self.x_min_traj = [];
                self.x_max_traj = [];
                self.t_min_traj = [];
                self.t_max_traj = [];
                self.v_meas_traj = [];
                self.r_meas_traj = [];
            else
                self.x_min_traj = Traj_con.x_min_traj;
                self.x_max_traj = Traj_con.x_max_traj;
                self.t_min_traj = Traj_con.t_min_traj;
                self.t_max_traj = Traj_con.r_max_traj;
                self.v_meas_traj = Traj_con.v_meas_traj;
                self.r_meas_traj = Traj_con.r_meas_traj;
            end
            
            % Density condition
            if isempty(Dens_con) || isempty(fieldnames(Dens_con))
                self.x_min_dens = [];
                self.x_max_dens = [];
                self.t_dens = [];
                self.dens_meas = [];
            else
                self.x_min_dens = Dens_con.x_min_dens;
                self.x_max_dens = Dens_con.x_max_dens;
                self.t_dens = Dens_con.t_dens;
                self.dens_meas = Dens_con.dens_meas;
            end
            
            
            % set initial condition
            if ~iscolumn(Initial_con.X_grid_cum)
                Initial_con.X_grid_cum = Initial_con.X_grid_cum';
            end
            self.X_grid_cum = Initial_con.X_grid_cum;
            
            if (length(Initial_con.X_grid_cum)~=length(Initial_con.IC)+1 ...
                && ~isempty(Initial_con.IC))
                errors('Check ini_segments and rho_ini dimensions');
            end
            
            self.X_grid = self.X_grid_cum(2:end) - self.X_grid_cum(1:end-1);
            
            if (isempty(Initial_con.IC) || ~all(~isnan(Initial_con.IC)))
                % This is because we have to check each initial condition,
                % instead of using heauristic knowledge and check only a
                % few conditins.
                warning(['Warning: Not defining a complete set of initial'...
                    'value conditions could increase the computational time!\n'])
            end
            
            if ~iscolumn(Initial_con.IC)
                Initial_con.IC = Initial_con.IC';
            end
            self.rho_ini = Initial_con.IC;
            
            % Use heuristic knowledge and group initial conditions into
            % blocks. By theory, only the first or the last in each block
            % need to be checked.
            if (~isempty(Initial_con.IC))
                block_ini = double(Initial_con.IC > self.k_c); 
                block_ini(isnan(Initial_con.IC)) = NaN;
                self.indiced_rho_ini = self.groupSameElement(block_ini,inf);
            else
                self.indiced_rho_ini = [];
            end
            
                      
            % Compute how many bool variables are needed if internal or
            % density conditions exist
            [self.nb_min_traj, self.nb_max_traj, self.nb_min_dens, self.nb_max_dens] = self.getBinaryvar;
            
            % index for bools
            % trajectory min points bools
            if sum(self.nb_min_traj) ~= 0
                self.dv_link.bool_internal_min = [self.dv_max+1;...
                                                  self.dv_max+sum(self.nb_min_traj)];
                self.dv_max = self.dv_max + sum(self.nb_min_traj);
            end
            % trajectory max points bools
            if sum(self.nb_max_traj) ~= 0
                self.dv_link.bool_internal_max = [self.dv_max+1;...
                                                  self.dv_max+sum(self.nb_max_traj)];
                self.dv_max = self.dv_max + sum(self.nb_max_traj);
            end
            % density min points bools
            if sum(self.nb_min_dens) ~= 0
                self.dv_link.bool_density_min = [self.dv_max+1;...
                                                  self.dv_max+sum(self.nb_min_dens)];
                self.dv_max = self.dv_max + sum(self.nb_min_dens);
            end
            % density max points bools
            if sum(self.nb_max_dens) ~= 0
                self.dv_link.bool_density_max = [self.dv_max+1;...
                                                  self.dv_max+sum(self.nb_max_dens)];
                self.dv_max = self.dv_max + sum(self.nb_max_dens);
            end
            
            % keep indexing all the bools 
            self.num_bool = sum(self.nb_min_traj) + sum(self.nb_max_traj) + ...
                            sum(self.nb_min_dens) + sum(self.nb_max_dens);
            if self.num_bool ~= 0
                self.dv_link.bool_all = [self.dv_var_max+1;...
                                         self.dv_var_max+self.num_bool];
            end
            
            % the last 1 is the right hand side constant for the
            % construction of the linear constraints
            self.size_row = self.dv_max + 1; 

        end
        
        
        
        %==========================================================================
        % Equations that define the upstream condition
        % input:
        %       t,x: the time, location of the point
        % output:
        %       array: an array which compute the value of the vehicle id
        %           using the upstream boundary conditions
        function[array] = us_con(self,t,x)
            
            array = zeros(1,self.size_row);  %Initialize the array
            
            if t < self.start_time || t > self.end_time
                array = [];
                return
            end
            
            n = sum(self.T_us_cum <= t);  % in step n interval [ , )
            if t == self.end_time     % the final time point
                n = n-1;    % in the last step interval
            end
            
            
            if( (abs(x - self.us_pos_m)<0.00001) && (self.T_us_cum(n)<=t) &&...
                (t<=self.T_us_cum(n+1)) && (n <= self.num_us_con) )
                array(1, self.dv_link.upstream(1) : self.dv_link.upstream(1)-1 + n-1) = self.T_us(1:n-1);
                array(1, self.dv_link.upstream(1)-1 + n) = t-self.T_us_cum(n);
                return
            else
                array = []; %return a "NULL"
                return
            end
        end
        
        
        %==========================================================================
        % Equations that define the downstream condition
        % input:
        %       t,x: the time, location of the point
        % output:
        %       array: an array which compute the value of the vehicle id
        %           using the downstream boundary conditions
        function [array] = ds_con(self,t,x)
            array = zeros(1,self.size_row);
            
            if t < self.start_time || t > self.end_time
                array = [];
                return
            end
            
            n = sum(self.T_ds_cum <= t);  % in step n interval [ , )
            if t == self.end_time     % the final time point
                n = n-1;    % in the last step interval
            end
            
            if( (abs(x-self.ds_pos_m)<0.01) && (self.T_ds_cum(n) <= t) &&...
                (t<= self.T_ds_cum(n+1)) && (n <= self.num_ds_con) )
                array(1,self.dv_link.initial(1):...
                        self.dv_link.initial(2) ) = -self.X_grid;   %Initial number of vehicles
                array(1,self.dv_link.downstream(1):self.dv_link.downstream(1)-1 + n-1) = self.T_ds(1:n-1);
                array(1,self.dv_link.downstream(1)-1 + n ) = t-self.T_ds_cum(n);
                return
            else
                array = [];
                return
            end
        end
        
        
        %==========================================================================
        % Equations that define the initial condition
        % input:
        %       b: the index of the initial condition this point in on.
        %       t,x: the time, location of the point
        % output:
        %       array: an array which compute the value of the vehicle id
        %           using the initial condition
        function[array] = initial_con(self, b, t, x)
            if( (abs(t-self.t0)<0.01) &&  ( sum(self.X_grid(1:b)) <= x) &&...
                    ( x <=sum(self.X_grid(1:b+1))) && (b <= self.num_initial_con) )
                array = zeros(1,self.size_row);
                array(1, self.dv_link.initial(1) : self.dv_link.initial(1)-1 + b) = -self.X_grid(1:b);
                array(1, self.dv_link.initial(1)-1 +  b +1) = -(x-sum(self.X_grid(1:b)));
                return
            else
                array = []; %Return a "NULL"
                return
            end
        end
        
        
        %==========================================================================
        % Equations that define the internal (trajectory) condition
        % input:
        %       m: the index of the internal condition this point in on.
        %       t,x: the time, location of the point
        % output:
        %       array: an array which compute the value of the vehicle id
        %           using the internal condition
        function [array] = traj_con(self, m ,t ,x)
            
            array = zeros(1,self.size_row);
            if ( (abs(x - (self.v_meas_traj(m)*(t-self.t_min_traj(m)) + self.x_min_traj(m))) < 0.01) &&...
                    (t>= self.t_min_traj(m)) && (t<= self.t_max_traj(m)) )
                array(1, self.dv_link.internal_L(1)-1 + m) = 1;
                array(1, self.dv_link.internal_r(1)-1 + m ) = t - self.t_min_traj(m);
                return
            else
                array = [];
                return
            end
        end
        
        
        %==========================================================================
        % Equations that define the density condition
        % input:
        %       u: the index of the density condition this point in on.
        %       t,x: the time, location of the point
        % output:
        %       array: an array which compute the value of the vehicle id
        %           using the density condition
        % new variables x_min_u, x_max_u, t_u, u_max, m_max_vm0
        function [array] = dens_con(self, u ,t ,x)
            
            if ( (abs(t-self.t_dens(u)) < 0.0001) &&...
                    (x>= self.x_min_dens(u)) && (x<= self.x_max_dens(u)) )
                array = zeros(1,self.size_row);
                array(1, self.dv_link.density_L(1)-1 + u) = 1;
                array(1, self.dv_link.density_rho(1)-1 + u ) =...
                    -(x - self.x_min_dens(u+1));
                return
            else
                
                array = [];
                return
                
            end
            
        end
        
        
        %==========================================================================
        % Equation that gives explicit vehicle id at (t,x) by the nth upstream condition
        % input:
        %       n: the index of upstream condition
        %       t,x: the time, location of the point
        % output:
        %       array: an array which compute the vehicle id at (t,x) using
        %           the nth upstream boundary condition
        function[array] = m_us_con(self,n,t,x)
            
            array = zeros(1,self.size_row);  %Initialize the array
            
            % characteristic domain
            if((self.T_us_cum(n) + (x - self.us_pos_m)/self.v <= t) && ...
                    (self.T_us_cum(n+1) + (x-self.us_pos_m)/self.v >= t ) && (n<= self.num_us_con))
                array(1, self.dv_link.upstream(1): self.dv_link.upstream(1)-1 + n-1) = self.T_us(1:n-1);
                array(1, self.dv_link.upstream(1)-1 + n) = t - (x-self.us_pos_m) / self.v - self.T_us_cum(n);
                return
                
            % Fan domain
            elseif((self.T_us_cum(n+1) + (x-self.us_pos_m)/self.v < t) && (n <= self.num_us_con))
                array(1,self.dv_link.upstream(1): self.dv_link.upstream(1)-1 + n) = self.T_us(1:n);
                array(1,self.size_row) = -self.k_c*self.v * (t - self.T_us_cum(n+1) - (x-self.us_pos_m)/self.v);
                return
                
            else
                array = [];     %return a "null"
                return
            end
            
        end
        
        
        %==========================================================================
        % Equation that gives explicit vehicle id at (t,x) by the nth downstream condition
        % input:
        %       n: the index of downstream condition
        %       t,x: the time, location of the point
        % output:
        %       array: an array which compute the vehicle id at (t,x) using
        %           the nth downstream boundary condition
        function [array] = m_ds_con(self,n,t,x)
            array = zeros(1,self.size_row);
            
            % in characteristic domain
            if( (self.T_ds_cum(n) + (x - self.ds_pos_m)/self.w <= t) &&...
                    (self.T_ds_cum(n+1) + (x-self.ds_pos_m)/self.w >= t )&& (n <= self.num_ds_con))
                array(1,self.dv_link.initial(1):...
                        self.dv_link.initial(2)) = -self.X_grid; %initial number of vehicles
                array(1,self.dv_link.downstream(1):self.dv_link.downstream(1)-1 + n-1) = self.T_ds(1:n-1);
                array(1,self.dv_link.downstream(1)-1 + n) = t - (x-self.ds_pos_m)/self.w - self.T_ds_cum(n);
                array(1,self.size_row) = self.k_m*(x-self.ds_pos_m);
                return
                
            % in Fan domain
            elseif( (self.T_ds_cum(n+1) + (x-self.ds_pos_m)/self.w < t) && (n <= self.num_ds_con))
                array(1,self.dv_link.initial(1):...
                        self.dv_link.initial(2)) = -self.X_grid; %initial number of vehicles
                array(1,self.dv_link.downstream(1):self.dv_link.downstream(1)-1 + n) = self.T_ds(1:n);
                array(1,self.size_row) = -self.k_c*self.v * (t - self.T_ds_cum(n+1) - (x-self.ds_pos_m)/(self.v));
                return
            else
                array = []; %return a "null"
                return
            end
        end
        
        
        %==========================================================================
        % Equation that gives explicit vehicle id at (t,x) by the mth
        %   internal (trajectory) condition
        % input:
        %       m: the index of internal condition
        %       t,x: the time, location of the point
        % output:
        %       array: an array which compute the vehicle id at (t,x) using
        %           the nth internal boundary condition
        function [array] = m_traj_con(self, m, t, x)
            array = zeros(1,self.size_row);
            
            if( (x >= self.x_min_traj(m) + self.v_meas_traj(m)*(t-self.t_min_traj(m))) &&...
                ( x >= self.x_max_traj(m) + self.v*(t-self.t_max_traj(m))) &&...
                (x <= self.x_min_traj(m) + self.v*(t-self.t_min_traj(m))) )
                array(1,self.dv_link.internal_L(1)-1 + m) = 1;
                array(1,self.dv_link.internal_r(1)-1 + m) =...
                    t - (x - self.x_min_traj(m) - self.v_meas_traj(m)*(t - self.t_min_traj(m))) / (self.v - self.v_meas_traj(m)) - self.t_min_traj(m) ;
                return
            end
            
            if ( (x<= self.x_min_traj(m) + self.v_meas_traj(m)*(t-self.t_min_traj(m))) &&...
                 ( x <= self.x_max_traj(m) + self.w*(t-self.t_max_traj(m))) &&...
                 (x >= self.x_min_traj(m) + self.w*(t-self.t_min_traj(m))) )
                array(1,self.dv_link.internal_L(1)-1 + m) = 1;
                array(1,self.dv_link.internal_r(1)-1 + m) =...
                    t - (x - self.x_min_traj(m) - self.v_meas_traj(m)*(t-self.t_min_traj(m)))/(self.w - self.v_meas_traj(m)) - self.t_min_traj(m);
                array(1,self.size_row) = -self.k_c*(self.v-self.w)*(x - self.x_min_traj(m) - self.v_meas_traj(m)*(t-self.t_min_traj(m)))/(self.w - self.v_meas_traj(m)) ;
                return
            end
            
            if( (x < self.x_max_traj(m) + self.v*(t - self.t_max_traj(m))) &&...
                    ( x > self.x_max_traj(m) + self.w*(t - self.t_max_traj(m))) )
                array(1, self.dv_link.internal_L(1)-1 + m) = 1;
                array(1, self.dv_link.internal_r(1)-1 + m) = self.t_max_traj(m) - self.t_min_traj(m);
                array(1, self.size_row) = -(t-self.t_max_traj(m))*self.k_c*(self.v - (x-self.x_max_traj(m))/(t-self.t_max_traj(m)+0.000001));
                return
            else
                array = []; %Return a "NULL"
                return
            end
        end
        
        
        %==========================================================================
        % Equation that gives explicit vehicle id at (t,x) by the bth
        %   initial condition, when the p<pc
        % input:
        %       b: the index of internal condition
        %       t,x: the time, location of the point
        % output:
        %       array: an array which compute the vehicle id at (t,x) using
        %           the nth internal boundary condition
        function [array] = m_initial_con_ff(self,b,t,x)
            array = zeros(1,self.size_row);
            
            % characteristic domain
            if ( (self.us_pos_m + sum(self.X_grid(1:b-1)) + t*self.v <= x) &&...
                    (x <=self.us_pos_m+ sum(self.X_grid(1:b)) + t*self.v) &&...
                    ( b <= self.num_initial_con))
                array(1, self.dv_link.initial(1) : ...
                         self.dv_link.initial(1)-1 + b-1) = -self.X_grid(1:b-1);
                
                array(1,self.dv_link.initial(1)-1 + b) = ...
                    (t*self.v - x + sum(self.X_grid(1:b-1)) + self.us_pos_m);
                return
            end
            
            % Fan domain
            if ( (self.us_pos_m+ sum(self.X_grid(1:b-1)) + t*self.w <= x) &&...
                    (x < self.us_pos_m+ sum(self.X_grid(1:b-1)) + t*self.v) &&...
                    (b <= self.num_initial_con))
                array(1, self.dv_link.initial(1): ...
                         self.dv_link.initial(1)-1 + b-1) = -self.X_grid(1:b-1);
                
                array(1,self.size_row) = -self.k_c*(t*self.v - x + sum(self.X_grid(1:b-1)) + self.us_pos_m);
                return
            else
                array = []; %Return a "NULL"
                return
            end
        end
        
        
        %==========================================================================
        % Equation that gives explicit vehicle id at (t,x) by the bth
        %   initial condition, when the p>pc
        % input:
        %       b: the index of initial condition
        %       t,x: the time, location of the point
        % output:
        %       array: an array which compute the vehicle id at (t,x) using
        %           the nth initial boundary condition
        function [array] = m_initial_con_cf(self,b,t,x)
            array = zeros(1,self.size_row);
            % characteristic domain
            if ( (self.us_pos_m + sum(self.X_grid(1:b-1)) + t*self.w <= x) &&...
                    (x <= self.us_pos_m + sum(self.X_grid(1:b)) + t*self.w) &&...
                    ( b <= self.num_initial_con))
                
                array(1,self.dv_link.initial(1) :...
                        self.dv_link.initial(1)-1 + b - 1 ) = -self.X_grid(1:b-1);
                array(1,self.dv_link.initial(1)-1 + b) = ...
                    (t*self.w - x + sum(self.X_grid(1:b-1)) + self.us_pos_m);
                array(1,self.size_row) = self.k_m*t*self.w;
                return
            end
            
            % Fan domain
            if ( (self.us_pos_m + sum(self.X_grid(1:b)) + t*self.w < x) &&...
                    (x <= self.us_pos_m + sum(self.X_grid(1:b)) + t*self.v) &&...
                    (b <= self.num_initial_con))
                
                array(1,self.dv_link.initial(1):...
                    self.dv_link.initial(1)-1 + b-1) = -self.X_grid(1:b-1);
                array(1,self.dv_link.initial(1)-1 + b) = -self.X_grid(b);
                array(1,self.size_row) = -(self.k_c*(t*self.w - x + sum(self.X_grid(1:b)) + self.us_pos_m) - self.k_m*t*self.w);
                return
            else
                array = []; %Return a "NULL" value
                return
            end
        end
        
        
        %==========================================================================        
        % Equation that gives explicit vehicle id at (t,x) by the uth
        %   density condition, when the p<pc
        % input:
        %       u: the index of density condition
        %       t,x: the time, location of the point
        % output:
        %       array: an array which compute the vehicle id at (t,x) using
        %           the nth density boundary condition
        function [array] = m_dens_con_ff(self,u,t,x)
            
            array = zeros(1,self.size_row);
            
            if ( (self.x_min_dens(u) + (t-self.t_dens(u))*self.v <= x) &&...
                 (x <=self.x_max_dens(u) + (t-self.t_dens(u))*self.v) &&...
                 (t>=self.t_dens(u)) && ( u <= self.num_density_con))
                
                array(1,self.dv_link.density_L(1)-1 + u ) = 1;
                array(1,self.dv_link.density_rho(1)-1 + u) =...
                    (t-self.t_dens(u))*self.v - x + self.x_min_dens(u);
                return
            end
            
            if ( (self.x_min_dens(u) + (t-self.t_dens(u))*self.w <= x) &&...
                 (x <= self.x_min_dens(u) + (t-self.t_dens(u))*self.v) &&...
                 (t>=self.t_dens(u)) && (u <= self.num_density_con))
                array(1,self.dv_link.density_L(1)-1 + u ) = 1;
                array(1,self.size_row) = -self.k_c*((t-self.t_dens(u))*self.v - x + self.x_min_dens(u));
                return
            else
                
                array = []; %Return a "NULL"
                return
                
            end
            
        end
      
        
        %==========================================================================        
        % Equation that gives explicit vehicle id at (t,x) by the uth
        %   density condition, when the p>pc
        % input:
        %       u: the index of density condition
        %       t,x: the time, location of the point
        % output:
        %       array: an array which compute the vehicle id at (t,x) using
        %           the nth density boundary condition
        function [array] = m_dens_con_cf(self,u,t,x)
            
            array = zeros(1,self.size_row);
            
            if ( (self.x_min_dens(u) + (t-self.t_dens(u))*self.w <= x) &&...
                 (x <=self.x_max_dens(u) + (t-self.t_dens(u))*self.w) &&...
                 (t>=self.t_dens(u)) && ( u <= self.num_density_con))
                
                array(1,self.dv_link.density_L(1)-1 + u ) = 1;
                array(1,self.dv_link.density_rho(1)-1 + u) =...
                        (t-self.t_dens(u))*self.w - x + self.x_min_dens(u);
                array(1,self.size_row) = self.k_m*(t-self.t_dens(u))*self.w;
                return
            end
            
            if ( (self.x_max_dens(u) + (t-self.t_dens(u))*self.w <= x) &&...
                 (x <= self.x_max_dens(u) + (t-self.t_dens(u))*self.v) &&...
                 (t>=self.t_dens(u)) && (u <= self.num_density_con))
                
                array(1,self.dv_link.density_L(1)-1 + u ) = 1;
                array(1,self.size_row) = -(self.k_c*((t-self.t_dens(u))*self.w - x +...
                                self.x_max_dens(u)) - self.k_m*(t-self.t_dens(u))*self.w);
                return
            else
                array = []; %Return a "NULL"
                return
            end
            
        end
        
        
        %==========================================================================
        % Function to create the model constraints matrix
        % Idea: the label at each point must be the smallest solution
        % output: 
        %       list[:,1:end-1] * x >= list[:,end]
        %       the model constraints matrix
        function  [list] = setModelMatrix(self)
            
            % first allocate memory for list
            
            list = zeros(100000,self.size_row);
            % initialize rows counts
            rows = 0; 
            
            %==============================================
            % Solution associated with n-th upstream boundary conditions
            for n=1:self.num_us_con
                
                % at downstream points including the begining and ending
                for p=1:self.num_ds_con+1
                    
                    if self.T_us_cum(n) + (self.ds_pos_m-self.us_pos_m)/self.v <= self.T_ds_cum(p) &&...
                       self.T_us_cum(n+1) + (self.ds_pos_m-self.us_pos_m)/self.v >= self.T_ds_cum(p)
                        % point (self.T_ds_cum(p), self.ds_pos_m)
                        array = self.m_us_con(n,self.T_ds_cum(p),self.ds_pos_m);
                        array2 = self.subtractArray(array, self.ds_con(self.T_ds_cum(p), self.ds_pos_m));
                        if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                            rows = rows+1;
                            list(rows,:) = array2;
                        end
                    end
                    
                end
                
                % freeflow speed intersection at downstream point
                % point (self.T_us_cum(n) + (self.ds_pos_m-self.us_pos_m)/self.v, self.ds_pos_m) 
                % where solution >= value condition 
                array = self.m_us_con(n, self.T_us_cum(n) + (self.ds_pos_m-self.us_pos_m)/self.v, self.ds_pos_m);
                array2 = self.subtractArray(array, self.ds_con(self.T_us_cum(n) + (self.ds_pos_m-self.us_pos_m)/self.v, self.ds_pos_m));
                if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                    rows = rows+1;
                    list(rows,:) = array2;
                end
                
                % There is no need to pose the inequality for upstream
                % points or initial time points
                
                % internal condition points
                for m=1:self.num_internal_con
                    
                    array = self.m_us_con(n,self.t_min_traj(m), self.x_min_traj(m));
                    array2 = self.subtractArray(array, self.traj_con(m, self.t_min_traj(m), self.x_min_traj(m)));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        
                        list(rows,:) = array2;
                    end
                    
                    array = self.m_us_con(n,self.t_max_traj(m), self.x_max_traj(m));
                    array2 = self.subtractArray(array, self.traj_con(m, self.t_max_traj(m), self.x_max_traj(m)));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.T_us_cum(n) * self.v - self.v_meas_traj(m) * self.t_min_traj(m) +...
                              self.x_min_traj(m) - self.us_pos_m) / (self.v - self.v_meas_traj(m));
                    x_temp = self.x_min_traj(m) + self.v_meas_traj(m)*(t_temp - self.t_min_traj(m));
                    array = self.m_us_con(n,t_temp, x_temp);
                    array2 = self.subtractArray(array, self.traj_con(m, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
          
                end
                
                % density condition points
                for u=1:self.num_density_con
                    
                    array = self.m_us_con(n,self.t_dens(u), self.x_min_dens(u));
                    array2 = self.subtractArray(array, self.dens_con(u, self.t_dens(u), self.x_min_dens(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    
                    array = self.m_us_con(n,self.t_dens(u), self.x_max_dens(u));
                    array2 = self.subtractArray(array, self.dens_con(u, self.t_dens(u), self.x_max_dens(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.m_us_con(n,self.t_dens(u), self.us_pos_m + self.v*(self.t_dens(u)-self.T_us_cum(n)));
                    array2 = self.subtractArray(array, self.dens_con(u, self.t_dens(u), self.us_pos_m +...
                        self.v*(self.t_dens(u)-self.T_us_cum(n))));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
     
                end
                
            end
            
            
            %==============================================
            % Solution associated with n-th downstram boundary conditions
            for n=1:self.num_ds_con
                
                % at upstream points
                for p=1:self.num_us_con+1
                    
                    if self.T_ds_cum(n) + (self.us_pos_m-self.ds_pos_m)/self.w <= self.T_us_cum(p) &&...
                            self.T_ds_cum(n+1) + (self.us_pos_m-self.ds_pos_m)/self.w >= self.T_us_cum(p)
                        % point (self.T_us_cum(p), self.us_pos_m)
                        array = self.m_ds_con(n,self.T_us_cum(p), self.us_pos_m);
                        array2 = self.subtractArray(array, self.us_con(self.T_us_cum(p), self.us_pos_m));
                        if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                            rows = rows+1;
                            
                            list(rows,:) = array2;
                        end
                    end
                    
                end
                
                % intersection point at upstream (self.T_ds_cum(n) + (self.us_pos_m-self.ds_pos_m)/self.w,self.us_pos_m)
                array = self.m_ds_con(n,self.T_ds_cum(n) + (self.us_pos_m-self.ds_pos_m)/self.w, self.us_pos_m);
                array2 = self.subtractArray(array, self.us_con(self.T_ds_cum(n) + (self.us_pos_m-self.ds_pos_m)/self.w, self.us_pos_m));
                if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                    rows = rows+1;
                    
                    list(rows,:) = array2;
                end
                
                % There is no need to pose the inequality for upstream
                % points or initial time points
                
                % internal condition points
                for m=1:self.num_internal_con
                    array = self.m_ds_con(n,self.t_min_traj(m), self.x_min_traj(m));
                    array2 = self.subtractArray(array, self.traj_con(m, self.t_min_traj(m), self.x_min_traj(m)));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.m_ds_con(n,self.t_max_traj(m), self.x_max_traj(m));
                    array2 = self.subtractArray(array, self.traj_con(m, self.t_max_traj(m), self.x_max_traj(m)));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.T_ds_cum(n) * self.w - self.v_meas_traj(m) * self.t_min_traj(m) +...
                        self.x_min_traj(m) - self.ds_pos_m) / (self.w - self.v_meas_traj(m));
                    x_temp = self.x_min_traj(m) + self.v_meas_traj(m)*(t_temp - self.t_min_traj(m));
                    array = self.m_ds_con(n,t_temp, x_temp);
                    array2 = self.subtractArray(array, self.traj_con(m, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                end
                
                % density condition points
                for u=1:self.num_density_con
                    
                    array = self.m_ds_con(n,self.t_dens(u), self.x_min_dens(u));
                    array2 = self.subtractArray(array, self.dens_con(u, self.t_dens(u), self.x_min_dens(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.m_ds_con(n,self.t_dens(u), self.x_max_dens(u));
                    array2 = self.subtractArray(array, self.dens_con(u, self.t_dens(u), self.x_max_dens(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.m_ds_con(n,self.t_dens(u), self.ds_pos_m + self.w*(self.t_dens(u)-self.T_ds_cum(n)));
                    array2 = self.subtractArray(array, self.dens_con(u, self.t_dens(u), self.ds_pos_m +...
                        self.w*(self.t_dens(u)-self.T_ds_cum(n))));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                end
                
            end
            
            
            %==============================================
            % Solution associated with initial conditions
            for k=1:self.num_initial_con
                
                % at upstream points
                for p=1:self.num_us_con+1
                    
                    %TAU1 rho < rho_c
                    
                    % upstream time grid points
                    array = self.m_initial_con_ff(k,self.T_us_cum(p),self.us_pos_m);
                    array2 = self.subtractArray(array,self.us_con(self.T_us_cum(p),self.us_pos_m));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    
                    %TAU 2
                    array = self.m_initial_con_cf(k,self.T_us_cum(p),self.us_pos_m);
                    array2 = self.subtractArray(array,self.us_con(self.T_us_cum(p),self.us_pos_m));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                end
                
                % at upstream intersection points
                % tau 1
                array = self.m_initial_con_ff(k,self.start_time + (-sum(self.X_grid(1:k)))/self.w,self.us_pos_m);
                array2 = self.subtractArray(array,self.us_con(self.start_time +...
                    (-sum(self.X_grid(1:k)))/self.w,self.us_pos_m));
                if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                    rows = rows+1;
                    list(rows,:) = array2;
                end
                
                % tau 2
                array = self.m_initial_con_cf(k,self.start_time + (-sum(self.X_grid(1:k)))/self.w,self.us_pos_m);
                array2 = self.subtractArray(array,self.us_con(self.start_time +...
                    (-sum(self.X_grid(1:k)))/self.w,self.us_pos_m));
                if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                    rows = rows+1;
                    list(rows,:) = array2;
                end
                
                
                % at downstream points
                for p=1:self.num_ds_con+1
                   
                    % Tau1
                    array = self.m_initial_con_ff(k,self.T_ds_cum(p),self.ds_pos_m);
                    array2 = self.subtractArray(array,self.ds_con(self.T_ds_cum(p),self.ds_pos_m));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                   
                    %tau2
                    array = self.m_initial_con_cf(k,self.T_ds_cum(p),self.ds_pos_m);
                    array2 = self.subtractArray(array,self.ds_con(self.T_ds_cum(p),self.ds_pos_m));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                end
                
                % intersection point at downstream
                array = self.m_initial_con_ff(k,self.start_time + (self.ds_pos_m-(sum(self.X_grid(1:k-1))+self.us_pos_m))/self.v,self.ds_pos_m);
                array2 = self.subtractArray(array,self.ds_con(self.start_time +...
                    (self.ds_pos_m-(sum(self.X_grid(1:k-1))+self.us_pos_m))/self.v,self.ds_pos_m));
                if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                    rows = rows+1;
                    list(rows,:) = array2;
                end
                
                array = self.m_initial_con_cf(k,self.start_time + (self.ds_pos_m-(sum(self.X_grid(1:k-1))+self.us_pos_m))/self.v,self.ds_pos_m);
                array2 = self.subtractArray(array,self.ds_con(self.start_time +...
                    (self.ds_pos_m-(sum(self.X_grid(1:k-1))+self.us_pos_m))/self.v,self.ds_pos_m));
                if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                    rows = rows+1;
                    list(rows,:) = array2;
                end
                
                % for internal points
                for p=1:self.num_internal_con
                    
                    %TAU 1
                    array = self.m_initial_con_ff(k,self.t_min_traj(p), self.x_min_traj(p));
                    array2 = self.subtractArray(array, self.traj_con(p, self.t_min_traj(p), self.x_min_traj(p)));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.m_initial_con_ff(k,self.t_max_traj(p), self.x_max_traj(p));
                    array2 = self.subtractArray(array, self.traj_con(p, self.t_max_traj(p), self.x_max_traj(p)));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.x_min_traj(p) - (self.us_pos_m + sum(self.X_grid(1:k))) - self.v_meas_traj(p) * self.t_min_traj(p)) ...
                        / (self.v - self.v_meas_traj(p));
                    x_temp = self.v_meas_traj(p) * (t_temp - self.t_min_traj(p)) + self.x_min_traj(p);
                    array = self.m_initial_con_ff(k,t_temp, x_temp);
                    array2 = self.subtractArray(array, self.traj_con(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.x_min_traj(p) - (self.us_pos_m + sum(self.X_grid(1:k-1))) - self.v_meas_traj(p) * self.t_min_traj(p))...
                        / (self.w - self.v_meas_traj(p));
                    x_temp = self.v_meas_traj(p) * (t_temp - self.t_min_traj(p)) + self.x_min_traj(p);
                    array = self.m_initial_con_ff(k,t_temp, x_temp);
                    array2 = self.subtractArray(array, self.traj_con(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.x_min_traj(p) - (self.us_pos_m + sum(self.X_grid(1:k-1))) - self.v_meas_traj(p) * self.t_min_traj(p))...
                        / (self.v - self.v_meas_traj(p));
                    x_temp = self.v_meas_traj(p) * (t_temp - self.t_min_traj(p)) + self.x_min_traj(p);
                    array = self.m_initial_con_ff(k,t_temp, x_temp);
                    array2 = self.subtractArray(array, self.traj_con(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.x_min_traj(p) - (self.us_pos_m + sum(self.X_grid(1:k))) - self.v_meas_traj(p) * self.t_min_traj(p))...
                        / (self.w - self.v_meas_traj(p));
                    x_temp = self.v_meas_traj(p) * (t_temp - self.t_min_traj(p)) + self.x_min_traj(p);
                    array = self.m_initial_con_ff(k,t_temp, x_temp);
                    array2 = self.subtractArray(array, self.traj_con(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    
                    %TAU 2
                    array = self.m_initial_con_cf(k,self.t_min_traj(p), self.x_min_traj(p));
                    array2 = self.subtractArray(array, self.traj_con(p, self.t_min_traj(p), self.x_min_traj(p)));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.m_initial_con_cf(k,self.t_max_traj(p), self.x_max_traj(p));
                    array2 = self.subtractArray(array, self.traj_con(p, self.t_max_traj(p), self.x_max_traj(p)));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.x_min_traj(p) - (self.us_pos_m + sum(self.X_grid(1:k))) - self.v_meas_traj(p) * self.t_min_traj(p))...
                        / (self.v - self.v_meas_traj(p));
                    x_temp = self.v_meas_traj(p) * (t_temp - self.t_min_traj(p)) + self.x_min_traj(p);
                    array = self.m_initial_con_cf(k,t_temp, x_temp);
                    array2 = self.subtractArray(array, self.traj_con(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.x_min_traj(p) - (self.us_pos_m + sum(self.X_grid(1:k-1))) - self.v_meas_traj(p) * self.t_min_traj(p))...
                        / (self.w - self.v_meas_traj(p));
                    x_temp = self.v_meas_traj(p) * (t_temp - self.t_min_traj(p)) + self.x_min_traj(p);
                    array = self.m_initial_con_cf(k,t_temp, x_temp);
                    array2 = self.subtractArray(array, self.traj_con(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.x_min_traj(p) - (self.us_pos_m + sum(self.X_grid(1:k-1))) - self.v_meas_traj(p) * self.t_min_traj(p))...
                        / (self.v - self.v_meas_traj(p));
                    x_temp = self.v_meas_traj(p) * (t_temp - self.t_min_traj(p)) + self.x_min_traj(p);
                    array = self.m_initial_con_cf(k,t_temp, x_temp);
                    array2 = self.subtractArray(array, self.traj_con(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.x_min_traj(p) - (self.us_pos_m + sum(self.X_grid(1:k))) - self.v_meas_traj(p) * self.t_min_traj(p))...
                        / (self.w - self.v_meas_traj(p));
                    x_temp = self.v_meas_traj(p) * (t_temp - self.t_min_traj(p)) + self.x_min_traj(p);
                    array = self.m_initial_con_cf(k,t_temp, x_temp);
                    array2 = self.subtractArray(array, self.traj_con(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                end
                
                % for density points
                for u=1:self.num_density_con
                    
                    %TAU 1
                    array = self.m_initial_con_ff(k,self.t_dens(u), self.x_min_dens(u));
                    array2 = self.subtractArray(array, self.dens_con(u, self.t_dens(u), self.x_min_dens(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.m_initial_con_ff(k,self.t_dens(u), self.x_max_dens(u));
                    array2 = self.subtractArray(array, self.dens_con(u, self.t_dens(u), self.x_max_dens(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.m_initial_con_ff(k,self.t_dens(u), self.us_pos_m + sum(self.X_grid(1:k-1)) +...
                        (self.t_dens(u)-self.start_time)*self.v);
                    array2 = self.subtractArray(array, self.dens_con(u, self.t_dens(u), self.us_pos_m + sum(self.X_grid(1:k-1)) +...
                        (self.t_dens(u)-self.start_time)*self.v));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.m_initial_con_ff(k,self.t_dens(u), self.us_pos_m + sum(self.X_grid(1:k)) +...
                        (self.t_dens(u)-self.start_time)*self.v);
                    array2 = self.subtractArray(array, self.dens_con(u, self.t_dens(u), self.us_pos_m + sum(self.X_grid(1:k)) +...
                        (self.t_dens(u)-self.start_time)*self.v));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.m_initial_con_ff(k,self.t_dens(u), self.us_pos_m + sum(self.X_grid(1:k-1)) +...
                        (self.t_dens(u)-self.start_time)*self.w);
                    array2 = self.subtractArray(array, self.dens_con(u, self.t_dens(u), self.us_pos_m + sum(self.X_grid(1:k-1)) +...
                        (self.t_dens(u)-self.start_time)*self.w));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.m_initial_con_ff(k,self.t_dens(u), self.us_pos_m + sum(self.X_grid(1:k)) +...
                        (self.t_dens(u)-self.start_time)*self.w);
                    array2 = self.subtractArray(array, self.dens_con(u, self.t_dens(u), self.us_pos_m + sum(self.X_grid(1:k)) +...
                        (self.t_dens(u)-self.start_time)*self.w));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    
                    %TAU 2
                    array = self.m_initial_con_cf(k,self.t_dens(u), self.x_min_dens(u));
                    array2 = self.subtractArray(array, self.dens_con(u, self.t_dens(u), self.x_min_dens(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.m_initial_con_cf(k,self.t_dens(u), self.x_max_dens(u));
                    array2 = self.subtractArray(array, self.dens_con(u, self.t_dens(u), self.x_max_dens(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.m_initial_con_cf(k,self.t_dens(u), self.us_pos_m + sum(self.X_grid(1:k-1)) +...
                        (self.t_dens(u)-self.start_time)*self.v);
                    array2 = self.subtractArray(array, self.dens_con(u, self.t_dens(u), self.us_pos_m + sum(self.X_grid(1:k-1)) +...
                        (self.t_dens(u)-self.start_time)*self.v));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.m_initial_con_cf(k,self.t_dens(u), self.us_pos_m + sum(self.X_grid(1:k)) +...
                        (self.t_dens(u)-self.start_time)*self.v);
                    array2 = self.subtractArray(array, self.dens_con(u, self.t_dens(u), self.us_pos_m + sum(self.X_grid(1:k)) +...
                        (self.t_dens(u)-self.start_time)*self.v));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.m_initial_con_cf(k,self.t_dens(u), self.us_pos_m + sum(self.X_grid(1:k-1)) +...
                        (self.t_dens(u)-self.start_time)*self.w);
                    array2 = self.subtractArray(array, self.dens_con(u, self.t_dens(u), self.us_pos_m + sum(self.X_grid(1:k-1)) +...
                        (self.t_dens(u)-self.start_time)*self.w));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.m_initial_con_cf(k,self.t_dens(u), self.us_pos_m + sum(self.X_grid(1:k)) +...
                        (self.t_dens(u)-self.start_time)*self.w);
                    array2 = self.subtractArray(array, self.dens_con(u, self.t_dens(u), self.us_pos_m + sum(self.X_grid(1:k)) +...
                        (self.t_dens(u)-self.start_time)*self.w));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                end
                
            end
            
            
            %==============================================
            % Solution associated with internal conditions
            for m=1:self.num_internal_con
                
                % upstream points
                for p=1:self.num_us_con+1
                    
                    array = self.m_traj_con(m,self.T_us_cum(p), self.us_pos_m);
                    array2 = self.subtractArray(array, self.us_con(self.T_us_cum(p), self.us_pos_m));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.us_pos_m - self.x_min_traj(m) + self.w*self.t_min_traj(m)) / (self.w);
                    x_temp = self.us_pos_m;
                    array = self.m_traj_con(m,t_temp, x_temp);
                    array2 = self.subtractArray(array, self.us_con(t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.us_pos_m - self.x_max_traj(m) + self.w*self.t_max_traj(m)) / (self.w);
                    x_temp = self.us_pos_m;
                    array = self.m_traj_con(m,t_temp, x_temp);
                    array2 = self.subtractArray(array, self.us_con(t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                end
                
                
                % downstream points
                for p=1:self.num_ds_con+1
                    
                    array = self.m_traj_con(m,self.T_ds_cum(p), self.ds_pos_m);
                    array2 = self.subtractArray(array, self.ds_con(self.T_ds_cum(p), self.ds_pos_m));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.ds_pos_m - self.x_min_traj(m) + self.v*self.t_min_traj(m)) / (self.v);
                    x_temp = self.ds_pos_m;
                    array = self.m_traj_con(m,t_temp, x_temp);
                    array2 = self.subtractArray(array, self.ds_con(t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.ds_pos_m - self.x_max_traj(m) + self.v*self.t_max_traj(m)) / (self.v);
                    x_temp = self.ds_pos_m;
                    array = self.m_traj_con(m,t_temp, x_temp);
                    array2 = self.subtractArray(array, self.ds_con(t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                end
                
                
                % other internal points
                for p=1:self.num_internal_con
                    
                    array = self.m_traj_con(m,self.t_min_traj(p), self.x_min_traj(p));
                    array2 = self.subtractArray(array, self.traj_con(p,self.t_min_traj(p), self.x_min_traj(p)));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.m_traj_con(m,self.t_max_traj(p), self.x_max_traj(p));
                    array2 = self.subtractArray(array, self.traj_con(p,self.t_max_traj(p), self.x_max_traj(p)));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.x_min_traj(m) - self.x_min_traj(p) + self.v_meas_traj(p) * self.t_min_traj(p) -...
                        self.v_meas_traj(m) * self.t_min_traj(m)) / (self.v_meas_traj(p) - self.v_meas_traj(m));
                    x_temp = self.v_meas_traj(p) * (t_temp - self.t_min_traj(p)) + self.x_min_traj(p);
                    array = self.m_traj_con(m,t_temp, x_temp);
                    array2 = self.subtractArray(array, self.traj_con(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    
                    t_temp = (self.x_max_traj(m) - self.x_min_traj(p) + self.v_meas_traj(p) * self.t_min_traj(p) -...
                        self.v * self.t_max_traj(m)) / (self.v_meas_traj(p) - self.v);
                    x_temp = self.v_meas_traj(p) * (t_temp - self.t_min_traj(p)) + self.x_min_traj(p);
                    array = self.m_traj_con(m,t_temp, x_temp);
                    array2 = self.subtractArray(array, self.traj_con(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    
                    t_temp = (self.x_min_traj(m) - self.x_min_traj(p) + self.v_meas_traj(p) * self.t_min_traj(p) -...
                        self.v * self.t_min_traj(m)) / (self.v_meas_traj(p) - self.v);
                    x_temp = self.v_meas_traj(p) * (t_temp - self.t_min_traj(p)) + self.x_min_traj(p);
                    array = self.m_traj_con(m,t_temp, x_temp);
                    array2 = self.subtractArray(array, self.traj_con(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    
                    t_temp = (self.x_max_traj(m) - self.x_min_traj(p) + self.v_meas_traj(p) * self.t_min_traj(p) -...
                        self.v * self.t_max_traj(m)) / (self.v_meas_traj(p) - self.w);
                    x_temp = self.v_meas_traj(p) * (t_temp - self.t_min_traj(p)) + self.x_min_traj(p);
                    array = self.m_traj_con(m,t_temp, x_temp);
                    array2 = self.subtractArray(array, self.traj_con(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    
                    t_temp = (self.x_min_traj(m) - self.x_min_traj(p) + self.v_meas_traj(p) * self.t_min_traj(p) -...
                        self.v * self.t_min_traj(m)) / (self.v_meas_traj(p) - self.w);
                    x_temp = self.v_meas_traj(p) * (t_temp - self.t_min_traj(p)) + self.x_min_traj(p);
                    array = self.m_traj_con(m,t_temp, x_temp);
                    array2 = self.subtractArray(array, self.traj_con(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
      
                    
                end
                
                
                % density condition points
                for u=1:self.num_density_con
                    
                    array = self.m_traj_con(m,self.t_dens(u), self.x_min_dens(u));
                    array2 = self.subtractArray(array, self.dens_con(u,self.t_dens(u), self.x_min_dens(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.m_traj_con(m,self.t_dens(u), self.x_max_dens(u));
                    array2 = self.subtractArray(array, self.dens_con(u,self.t_dens(u), self.x_max_dens(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.m_traj_con(m,self.t_dens(u), self.x_min_traj(m)+(self.t_dens(u)-self.t_min_traj(m))*self.w);
                    array2 = self.subtractArray(array, self.dens_con(u,self.t_dens(u),...
                                               self.x_min_traj(m)+(self.t_dens(u)-self.t_min_traj(m))*self.w));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.m_traj_con(m,self.t_dens(u), self.x_max_traj(m)+(self.t_dens(u)-self.t_max_traj(m))*self.w);
                    array2 = self.subtractArray(array, self.dens_con(u,self.t_dens(u),...
                        self.x_max_traj(m)+(self.t_dens(u)-self.t_max_traj(m))*self.w));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.m_traj_con(m,self.t_dens(u), self.x_min_traj(m)+(self.t_dens(u)-self.t_min_traj(m))*self.v);
                    array2 = self.subtractArray(array, self.dens_con(u,self.t_dens(u),...
                        self.x_min_traj(m)+(self.t_dens(u)-self.t_min_traj(m))*self.v));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.m_traj_con(m,self.t_dens(u), self.x_max_traj(m)+(self.t_dens(u)-self.t_max_traj(m))*self.v);
                    array2 = self.subtractArray(array, self.dens_con(u,self.t_dens(u),...
                        self.x_max_traj(m)+(self.t_dens(u)-self.t_max_traj(m))*self.v));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                end
                
            end
            
            
            % Solution associated with density conditions
            for k=1:self.num_density_con

                % upstream points
                for p=1:self.num_us_con+1
                    
                    %UPS1
                    array = self.m_dens_con_ff(k,self.T_us_cum(p),self.us_pos_m);
                    array2 = self.subtractArray(array,self.us_con(self.T_us_cum(p),self.us_pos_m));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    %UPS2
                    array = self.m_dens_con_cf(k,self.T_us_cum(p),self.us_pos_m);
                    array2 = self.subtractArray(array,self.us_con(self.T_us_cum(p),self.us_pos_m));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                end
                
                % upstream intersection points
                array = self.m_dens_con_ff(k,self.t_dens(k) + (self.us_pos_m-self.x_max_dens(k))/self.w,self.us_pos_m);
                array2 = self.subtractArray(array,self.us_con(self.t_dens(k) +...
                        (self.us_pos_m-self.x_max_dens(k))/self.w,self.us_pos_m));
                if(~isempty(array2))
                    rows = rows+1;
                    list(rows,:) = array2;
                end
                
                array = self.m_dens_con_cf(k,self.t_dens(k) + (self.us_pos_m-self.x_max_dens(k))/self.w,self.us_pos_m);
                array2 = self.subtractArray(array,self.us_con(self.t_dens(k) +...
                    (self.us_pos_m-self.x_max_dens(k))/self.w,self.us_pos_m));
                if(~isempty(array2))
                    rows = rows+1;
                    list(rows,:) = array2;
                end
                
                
                % downstream points
                for p=1:self.num_ds_con+1
                    
                    %UPS1
                    array = self.m_dens_con_ff(k,self.T_ds_cum(p),self.ds_pos_m);
                    array2 = self.subtractArray(array,self.ds_con(self.T_ds_cum(p),self.ds_pos_m));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    %USP2
                    array = self.m_dens_con_cf(k,self.T_ds_cum(p),self.ds_pos_m);
                    array2 = self.subtractArray(array,self.ds_con(self.T_ds_cum(p),self.ds_pos_m));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                end
                
                % downstream intersection points
                array = self.m_dens_con_ff(k,self.t_dens(k) + (self.ds_pos_m-self.x_min_dens(k))/self.v,self.ds_pos_m);
                array2 = self.subtractArray(array,self.ds_con(self.t_dens(k) +...
                    (self.ds_pos_m-self.x_min_dens(k))/self.v,self.ds_pos_m));
                if(~isempty(array2))
                    rows = rows+1;
                    list(rows,:) = array2;
                end
                
                array = self.m_dens_con_cf(k,self.t_dens(k) + (self.ds_pos_m-self.x_min_dens(k))/self.v,self.ds_pos_m);
                array2 = self.subtractArray(array,self.ds_con(self.t_dens(k) +...
                    (self.ds_pos_m-self.x_min_dens(k))/self.v,self.ds_pos_m));
                if(~isempty(array2))
                    rows = rows+1;
                    list(rows,:) = array2;
                end
                
                
                % internal condition points
                for p=1:self.num_internal_con
                    
                    %UPS 1
                    array = self.m_dens_con_ff(k,self.t_min_traj(p), self.x_min_traj(p));
                    array2 = self.subtractArray(array, self.traj_con(p, self.t_min_traj(p), self.x_min_traj(p)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.m_dens_con_ff(k,self.t_max_traj(p), self.x_max_traj(p));
                    array2 = self.subtractArray(array, self.traj_con(p, self.t_max_traj(p), self.x_max_traj(p)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.x_min_traj(p) - self.x_max_dens(k) - self.v_meas_traj(p) * self.t_min_traj(p) +...
                        self.v*self.t_dens(k)) / (self.v - self.v_meas_traj(p));
                    x_temp = self.v_meas_traj(p) * (t_temp - self.t_min_traj(p)) + self.x_min_traj(p);
                    array = self.m_dens_con_ff(k,t_temp, x_temp);
                    array2 = self.subtractArray(array, self.traj_con(p, t_temp, x_temp));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.x_min_traj(p) - self.x_min_dens(k) - self.v_meas_traj(p) * self.t_min_traj(p) +...
                        self.w*self.t_dens(k)) / (self.w - self.v_meas_traj(p));
                    x_temp = self.v_meas_traj(p) * (t_temp - self.t_min_traj(p)) + self.x_min_traj(p);
                    array = self.m_dens_con_ff(k,t_temp, x_temp);
                    array2 = self.subtractArray(array, self.traj_con(p, t_temp, x_temp));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.x_min_traj(p) - self.x_min_dens(k) - self.v_meas_traj(p) * self.t_min_traj(p) +...
                        self.v*self.t_dens(k)) / (self.v - self.v_meas_traj(p));
                    x_temp = self.v_meas_traj(p) * (t_temp - self.t_min_traj(p)) + self.x_min_traj(p);
                    array = self.m_dens_con_ff(k,t_temp, x_temp);
                    array2 = self.subtractArray(array, self.traj_con(p, t_temp, x_temp));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.x_min_traj(p) - self.x_max_dens(k) - self.v_meas_traj(p) * self.t_min_traj(p) +...
                        self.w*self.t_dens(k)) / (self.w - self.v_meas_traj(p));
                    x_temp = self.v_meas_traj(p) * (t_temp - self.t_min_traj(p)) + self.x_min_traj(p);
                    array = self.m_dens_con_ff(k,t_temp, x_temp);
                    array2 = self.subtractArray(array, self.traj_con(p, t_temp, x_temp));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    
                    %UPS 2
                    array = self.m_dens_con_cf(k,self.t_min_traj(p), self.x_min_traj(p));
                    array2 = self.subtractArray(array, self.traj_con(p, self.t_min_traj(p), self.x_min_traj(p)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.m_dens_con_cf(k,self.t_max_traj(p), self.x_max_traj(p));
                    array2 = self.subtractArray(array, self.traj_con(p, self.t_max_traj(p), self.x_max_traj(p)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.x_min_traj(p) - self.x_max_dens(k) - self.v_meas_traj(p) * self.t_min_traj(p) +...
                        self.v*self.t_dens(k)) / (self.v - self.v_meas_traj(p));
                    x_temp = self.v_meas_traj(p) * (t_temp - self.t_min_traj(p)) + self.x_min_traj(p);
                    array = self.m_dens_con_cf(k,t_temp, x_temp);
                    array2 = self.subtractArray(array, self.traj_con(p, t_temp, x_temp));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.x_min_traj(p) - self.x_min_dens(k) - self.v_meas_traj(p) * self.t_min_traj(p) +...
                        self.w*self.t_dens(k)) / (self.w - self.v_meas_traj(p));
                    x_temp = self.v_meas_traj(p) * (t_temp - self.t_min_traj(p)) + self.x_min_traj(p);
                    array = self.m_dens_con_cf(k,t_temp, x_temp);
                    array2 = self.subtractArray(array, self.traj_con(p, t_temp, x_temp));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.x_min_traj(p) - self.x_min_dens(k) - self.v_meas_traj(p) * self.t_min_traj(p) +...
                        self.v*self.t_dens(k)) / (self.v - self.v_meas_traj(p));
                    x_temp = self.v_meas_traj(p) * (t_temp - self.t_min_traj(p)) + self.x_min_traj(p);
                    array = self.m_dens_con_cf(k,t_temp, x_temp);
                    array2 = self.subtractArray(array, self.traj_con(p, t_temp, x_temp));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.x_min_traj(p) - self.x_max_dens(k) - self.v_meas_traj(p) * self.t_min_traj(p) +...
                        self.w*self.t_dens(k)) / (self.w - self.v_meas_traj(p));
                    x_temp = self.v_meas_traj(p) * (t_temp - self.t_min_traj(p)) + self.x_min_traj(p);
                    array = self.m_dens_con_cf(k,t_temp, x_temp);
                    array2 = self.subtractArray(array, self.traj_con(p, t_temp, x_temp));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                end
                
                
                % density condition points
                for u=1:self.num_density_con
                    
                    %UPS 1
                    array = self.m_dens_con_ff(k,self.t_dens(u), self.x_min_dens(u));
                    array2 = self.subtractArray(array, self.dens_con(u, self.t_dens(u), self.x_min_dens(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.m_dens_con_ff(k,self.t_dens(u), self.x_max_dens(u));
                    array2 = self.subtractArray(array, self.dens_con(u, self.t_dens(u), self.x_max_dens(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.m_dens_con_ff(k,self.t_dens(u), self.x_min_dens(k) + (self.t_dens(u)-self.t_dens(k))*self.v);
                    array2 = self.subtractArray(array, self.dens_con(u, self.t_dens(u), self.x_min_dens(k) +...
                        (self.t_dens(u)-self.t_dens(k))*self.v));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.m_dens_con_ff(k,self.t_dens(u), self.x_max_dens(k) + (self.t_dens(u)-self.t_dens(k))*self.v);
                    array2 = self.subtractArray(array, self.dens_con(u, self.t_dens(u), self.x_max_dens(k) +...
                        (self.t_dens(u)-self.t_dens(k))*self.v));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.m_dens_con_ff(k,self.t_dens(u), self.x_min_dens(k) + (self.t_dens(u)-self.t_dens(k))*self.w);
                    array2 = self.subtractArray(array, self.dens_con(u, self.t_dens(u), self.x_min_dens(k) +...
                        (self.t_dens(u)-self.t_dens(k))*self.w));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.m_dens_con_ff(k,self.t_dens(u), self.x_max_dens(k) + (self.t_dens(u)-self.t_dens(k))*self.w);
                    array2 = self.subtractArray(array, self.dens_con(u, self.t_dens(u), self.x_max_dens(k) +...
                        (self.t_dens(u)-self.t_dens(k))*self.w));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    % %UPS 2
                    array = self.m_dens_con_cf(k,self.t_dens(u), self.x_min_dens(u));
                    array2 = self.subtractArray(array, self.dens_con(u, self.t_dens(u), self.x_min_dens(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.m_dens_con_cf(k,self.t_dens(u), self.x_max_dens(u));
                    array2 = self.subtractArray(array, self.dens_con(u, self.t_dens(u), self.x_max_dens(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.m_dens_con_cf(k,self.t_dens(u), self.x_min_dens(k) + (self.t_dens(u)-self.t_dens(k))*self.v);
                    array2 = self.subtractArray(array, self.dens_con(u, self.t_dens(u), self.x_min_dens(k) +...
                        (self.t_dens(u)-self.t_dens(k))*self.v));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.m_dens_con_cf(k,self.t_dens(u), self.x_max_dens(k) + (self.t_dens(u)-self.t_dens(k))*self.v);
                    array2 = self.subtractArray(array, self.dens_con(u, self.t_dens(u), self.x_max_dens(k) +...
                        (self.t_dens(u)-self.t_dens(k))*self.v));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.m_dens_con_cf(k,self.t_dens(u), self.x_min_dens(k) + (self.t_dens(u)-self.t_dens(k))*self.w);
                    array2 = self.subtractArray(array, self.dens_con(u, self.t_dens(u), self.x_min_dens(k) +...
                        (self.t_dens(u)-self.t_dens(k))*self.w));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.m_dens_con_cf(k,self.t_dens(u), self.x_max_dens(k) + (self.t_dens(u)-self.t_dens(k))*self.w);
                    array2 = self.subtractArray(array, self.dens_con(u, self.t_dens(u), self.x_max_dens(k) +...
                        (self.t_dens(u)-self.t_dens(k))*self.w));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                end
                
            end
            
            
            % truncate the list matrix to the first rows with
            % inequality information.
            list = list(1:rows,:);
            
            
        end
        
        
        %==========================================================================
        % Function to create the data constraints
        % Relative error does not make sense and fail when x_meas = 0.
        % Use absolute error expressed by percent of the normal range.
        % The corresponding variable x is constrained by
        % x - x_normal*e <= x <= x + x_normal*e, where x_normal = q_max, kc
        % output: 
        %       list[:,1:end-1] * x >= list[:,end]
        %       the data constraints matrix
        function [list2] = setDataMatrix(self)
                        
            %Use sparse matrix to speed up computation
            list2 = zeros(10000,self.size_row);
            rows = 0;   %initialize row counts
            
            % upstream data constraints
            if (~isempty(self.qin_meas))
                
                for n=1:self.num_us_con
                    if ~isnan(self.qin_meas(n))
                        array = zeros(1,self.size_row);
                        array(1,self.dv_link.upstream(1)-1 + n) = 1;
                        err = (self.T_us_cum(n+1)<= self.now_time)*self.e_meas_flow +...
                              (self.T_us_cum(n+1)> self.now_time)*self.e_his;
                        array(1,self.size_row) = self.qin_meas(n)-err*self.q_max;
                        
                        rows = rows+1;
                        list2(rows,:)=array;
                    end
                end
                
                for n=1:self.num_us_con
                    if ~isnan(self.qin_meas(n))
                        array = zeros(1,self.size_row);
                        array(1,self.dv_link.upstream(1)-1 + n) = -1;
                        err = (self.T_us_cum(n+1)<=self.now_time)*self.e_meas_flow +...
                            (self.T_us_cum(n+1)> self.now_time)*self.e_his;
                        array(1,self.size_row) = -self.qin_meas(n)-err*self.q_max;
                        rows = rows+1;
                        list2(rows,:)=array;
                    end
                end
                
            end
            
            % downstream data constraints
            if (~isempty(self.qout_meas) )
                
                for n=1:self.num_ds_con
                    if ~isnan(self.qout_meas(n))
                    array = zeros(1,self.size_row);
                    array(1,self.dv_link.downstream(1)-1 + n) = 1;
                    err = (self.T_ds_cum(n+1)<=self.now_time)*self.e_meas_flow +...
                              (self.T_ds_cum(n+1) > self.now_time)*self.e_his;
                    array(1,self.size_row) = self.qout_meas(n)-err*self.q_max;
                    rows = rows+1;
                    list2(rows,:)=array;
                    end
                end
                
                for n=1:self.num_ds_con
                    if ~isnan(self.qout_meas(n))
                    array = zeros(1,self.size_row);
                    array(1,self.dv_link.downstream(1)-1 + n) = -1;
                    err = (self.T_ds_cum(n+1)<=self.now_time)*self.e_meas_flow +...
                              (self.T_ds_cum(n+1)> self.now_time)*self.e_his;
                    array(1,self.size_row) = -self.qout_meas(n)-err*self.q_max;
                    rows = rows+1;
                    list2(rows,:)=array;
                    end
                end
                
            end
            
            % Initial data constraints
            if(~isempty(self.rho_ini) )
                
                for n=1:self.num_initial_con
                    if ~isnan(self.rho_ini(n))
                        array = zeros(1,self.size_row);
                        array(1,self.dv_link.initial(1)-1 + n) = 1;
                        array(1,self.size_row) = self.rho_ini(n) - self.e_est*self.k_c;
                        rows = rows+1;
                        list2(rows,:)=array;
                    end
                end
                
                for n=1:self.num_initial_con
                    if ~isnan(self.rho_ini(n))
                        array = zeros(1,self.size_row);
                        array(1,self.dv_link.initial(1)-1 + n) = -1;
                        array(1,self.size_row) = -self.rho_ini(n) - self.e_est*self.k_c;
                        rows = rows+1;
                        list2(rows,:)=array;
                    end
                end
            end
            
            % The relative lable of the vehicle may not be know, but we can
            % roughly know the passing rate, set here. Since no knowledge
            % on the normal range of r, use relative error
            % internal data constraints for rate of change r
            if(~isempty(self.r_meas_traj))
                
                for m = 1:self.num_internal_con
                    if ~isnan(self.r_meas_traj(m))
                        array = zeros(1,self.size_row);
                        array(1,self.dv_link.internal_r(1)-1 + m) = 1;
                        array(1,self.size_row) = self.r_meas_traj(m)*(1-self.e_meas_traj_r);
                        rows = rows+1;
                        list2(rows,:)=array;
                    end
                end
                
                for m = 1:self.num_internal_con
                    if ~isnan(self.r_meas_traj(m))
                        array = zeros(1,self.size_row);
                        array(1,self.dv_link.internal_r(1)-1 + m) = -1;
                        array(1,self.size_row) = - self.r_meas_traj(m)*(1+self.e_meas_traj_r);
                        rows = rows+1;
                        list2(rows,:)=array;
                    end
                end
                
            end
            
            % density data constraints for density rho_u
            if(~isempty(self.dens_meas))
                
                for u = 1:self.num_density_con
                    if ~isnan(self.dens_meas(u))
                        array = zeros(1,self.size_row);
                        array(1,self.dv_link.density_rho(1)-1 + u) = 1;
                        array(1,self.size_row) = self.dens_meas(u)-self.e_meas_dens*self.k_c;
                        rows = rows+1;
                        list2(rows,:)=array;
                    end
                end
                
                for u = 1:self.num_density_con
                    if ~isnan(self.dens_meas(u))
                        array = zeros(1,self.size_row);
                        array(1,self.dv_link.density_rho(1)-1 + u) = -1;
                        array(1,self.size_row) = -self.dens_meas(u)-self.e_meas_dens*self.k_c;
                        rows = rows+1;
                        list2(rows,:)=array;
                    end
                end
                
            end
            
            % truncate matrix
            list2 = list2(1:rows,:);
            
            
        end
        
        
        
        %==========================================================================
        % Function to create the MILP constraints associated with internal
        % trajectory conditions. Each internal trajectory point has a
        % vehicle id L.
        % L = min( Mi (all value condition solutions) )
        % Use binary variables to remove nonliear min function as
        %   L <= Mi for all i
        %   inf*bi + L >= Mi for all i
        %   sum(1-bi) = 1 (one and only one bi = 0)
        % output: 
        %       list[:,1:end-1] * x >= list[:,end]
        function [list3] = setInternalConstraints(self)
            
            list3 = zeros(0,0);
            
            %Binary variables useful counter as a reference
            countB=0;
            Cmax = inf; 
            
            % nb: number of binary variables per point
            nb = ceil(log2(2*(self.num_initial_con) + ...
                    (self.num_us_con)+ (self.num_ds_con) +...
                    (self.num_internal_con)+(self.num_density_con)));
            
            % comb_mat: combination matrix for possible binary combinations
            comb_mat = zeros(2^nb,nb);
            
            for i = 1:2^nb
                comb_mat(i,:) = double(dec2bin(i-1,nb)=='1');
            end
            
            % for each internal trajectory condtion, we set min point first
            % and then max point
            for k = 1:self.num_internal_con
                
                % ------------------------------------------------------------------------
                % x_min and t_min point------------------
                if (k==1 || (self.x_min_traj(k)~=self.x_max_traj(k-1) &&...
                        self.t_min_traj(k)~=self.t_max_traj(k-1)))
                    tempM = zeros(0,0);
                    %Upper bound constraints
                    %L1<=M1
                    %L1<=M2
                    % The internal trajectory condition value Lk
                    temp_array = self.traj_con(k, self.t_min_traj(k), self.x_min_traj(k));
                    
                    % Initial condition solutions Mi
                    % Lk <= Mi
                    for block = 1:size(self.indiced_rho_ini,1) 
                        
                        if (self.indiced_rho_ini(block,3) == 0)
                            
                            for b = self.indiced_rho_ini(block,1):self.indiced_rho_ini(block,2)
                                %-1 to be consistent with other functions self.mtau...
                                arraym = self.m_initial_con_ff(b,self.t_min_traj(k), self.x_min_traj(k));
                                array2 = self.subtractArray(arraym,temp_array);
                                if(~isempty(array2))
                                    rows = size(list3,1);
                                    list3(rows+1,:) = array2;
                                    rows = size(tempM,1);
                                    tempM(rows+1,:) = arraym;
                                    break;  % Only need to consider the first one
                                end
                            end
                            
                        elseif (self.indiced_rho_ini(block,3) == 1)
                            
                            for b = self.indiced_rho_ini(block,2):-1:self.indiced_rho_ini(block,1)
                                % For k>k_c case, From the last one
                                arraym = self.m_initial_con_cf(b,self.t_min_traj(k), self.x_min_traj(k));
                                array2 = self.subtractArray(arraym,temp_array);
                                if(~isempty(array2))
                                    rows = size(list3,1);
                                    list3(rows+1,:) = array2;
                                    rows = size(tempM,1);
                                    tempM(rows+1,:) = arraym;
                                    break;
                                end
                            end
                            
                        elseif (isnan(self.indiced_rho_ini(block,3)))
                            
                            for b = self.indiced_rho_ini(block,1):self.indiced_rho_ini(block,2)
                                % For undefined initial blocks, check all
                                arraym = self.m_initial_con_ff(b,self.t_min_traj(k), self.x_min_traj(k));
                                array2 = self.subtractArray(arraym,temp_array);
                                if(~isempty(array2))
                                    rows = size(list3,1);
                                    list3(rows+1,:) = array2;
                                    rows = size(tempM,1);
                                    tempM(rows+1,:) = arraym;
                                end
                                
                                arraym = self.m_initial_con_cf(b,self.t_min_traj(k), self.x_min_traj(k));
                                array2 = self.subtractArray(arraym,temp_array);
                                if(~isempty(array2))
                                    rows = size(list3,1);
                                    list3(rows+1,:) = array2;
                                    rows = size(tempM,1);
                                    tempM(rows+1,:) = arraym;
                                end
                            end
                            
                        end
                    end 
                    
                
                    
                    % Downstream condition solutions Mi
                    % Lk <= Mi
                    for n=1:self.num_ds_con
                        
                        arraym = self.m_ds_con(n,self.t_min_traj(k), self.x_min_traj(k));
                        array2 = self.subtractArray(arraym,temp_array);
                        if(~isempty(array2) && arraym(1, self.dv_link.downstream(1)-1 + n) ~= self.T_ds(n))
                            % The second condition ensure we only consider
                            % the solution at the characteristic domain
                            %Add a constraint
                            rows = size(list3,1);
                            list3(rows+1,:) = array2;
                            %Assign the solution to temporary Matrix
                            rows = size(tempM,1);
                            tempM(rows+1,:) = arraym;
                            break;  %only one
                        end
                        
                    end
                    
                    % Upstream condition solutions Mi
                    % Lk <= Mi
                    for n=1:self.num_us_con    
                        
                        arraym = self.m_us_con(n,self.t_min_traj(k), self.x_min_traj(k));
                        array2 = self.subtractArray(arraym, temp_array);
                        if(~isempty(array2) && arraym(1,self.size_row)==0)
                            %Add a constraint
                            rows = size(list3,1);
                            list3(rows+1,:) = array2;
                            %Assign the solution to temporary Matrix
                            rows = size(tempM,1);
                            tempM(rows+1,:) = arraym;
                            break;
                        end
                        
                    end
                    
                    % Internal condition solutions Mi
                    % Lk <= Mi
                    for n=1:self.num_internal_con   
                        
                        if(n~=k)
                            arraym = self.m_traj_con(n,self.t_min_traj(k), self.x_min_traj(k));
                            array2 = self.subtractArray(arraym,temp_array);
                            if(~isempty(array2))
                                rows = size(list3,1);
                                list3(rows+1,:) = array2;
                                rows = size(tempM,1);
                                tempM(rows+1,:) = arraym;
                            end
                        end
                        
                    end
                    
                    % Density condition solutions Mi
                    % Lk <= Mi
                    % Find density conditions that are not at t_end. Not need to compute solutions
                    % associated with density conditions at t_end
                    index_u = find(self.t_dens ~= self.end_time);
                    for i=1:length(index_u) % For density
                        
                        % for each density condition
                        n = index_u(i);
                        
                        arraym = self.m_dens_con_ff(n,self.t_min_traj(k), self.x_min_traj(k));
                        array2 = self.subtractArray(arraym,temp_array);
                        if(~isempty(array2))
                            rows = size(list3,1);
                            list3(rows+1,:) = array2;
                            rows = size(tempM,1);
                            tempM(rows+1,:) = arraym;
                        end
                        
                        arraym = self.m_dens_con_cf(n,self.t_min_traj(k), self.x_min_traj(k));
                        array2 = self.subtractArray(arraym,temp_array);
                        if(~isempty(array2))
                            rows = size(list3,1);
                            list3(rows+1,:) = array2;
                            rows = size(tempM,1);
                            tempM(rows+1,:) = arraym;
                        end
                        
                    end
                    
                    rowsM = size(tempM,1);
                    

                    % Define all the possible combinations according to the solutions that apply
                    % at the specified point
                    % Lower bound constraints
                    % Each row of omb_mat select one solution, for instance (b1, b2) = (1, 0)
                    % (1-b1)*inf + b2*inf + L >= M2 selects M2 solution.                     
                    for i = 1:rowsM
                        
                        array = temp_array;
                        
                        %Decode the binary combinations
                        for counter=1:self.nb_min_traj(k)
                            if comb_mat(i,nb-self.nb_min_traj(k) + counter) == 1
                                array(1, self.dv_link.bool_internal_min(1)-1 + countB + counter) = -Cmax;
                                
                            elseif comb_mat(i,nb-self.nb_min_traj(k) + counter) == 0
                                array(1, self.dv_link.bool_internal_min(1)-1 + countB + counter) = Cmax;
                                
                            end
                        end
                        array(1,self.size_row) = -Cmax*(sum(comb_mat(i,:))); %RHS (set negative to be on same side)
                        
                        % Add constraint to the MILP matrix
                        array2 = self.subtractArray(array,tempM(i,:));
                        rows = size(list3,1);
                        list3(rows+1,:) = array2;
                    end
                    
                    % Define the last constraint to bound the binary combination
                    % We may have only 3 solutions, while we used 2 binary variabels to describe the combination. 
                    % Then we bound 2*b1+1*b2 <= 3 the limit the binary search in those meaningful combinations.
                    array = zeros(1,self.size_row);
                    
                    for counter = 1:self.nb_min_traj(k)
                        
                        array(1, self.dv_link.bool_internal_min(1)-1 + countB + counter) = -2^(self.nb_min_traj(k)-counter);
                        
                    end
                    array(1,self.size_row) = -(rowsM-1); %RHS (maximum possible value of the binary comb)                    
                    rows = size(list3,1);
                    list3(rows+1,:) = array;
                    
                    countB = countB+self.nb_min_traj(k);
                    
                end
                
                % ------------------------------------------------------------------------
                % x_max and t_max point-------------------------
                tempM = zeros(0,0);
                
                %Upper bound constraints
                %L1<=M1
                %L1<=M2
                temp_array = self.traj_con(k, self.t_max_traj(k), self.x_max_traj(k));
                
                % Solutions associated with initial conditions
                for block = 1:size(self.indiced_rho_ini,1) 
                    
                    if (self.indiced_rho_ini(block,3) == 0)
                        
                        for b = self.indiced_rho_ini(block,1):self.indiced_rho_ini(block,2)
                            %-1 to be consistent with other functions self.mtau...
                            arraym = self.m_initial_con_ff(b,self.t_max_traj(k), self.x_max_traj(k));
                            array2 = self.subtractArray(arraym,temp_array);
                            if(~isempty(array2))
                                rows = size(list3,1);
                                list3(rows+1,:) = array2;
                                rows = size(tempM,1);
                                tempM(rows+1,:) = arraym;
                                break;  % Only need to consider the first one
                            end
                        end
                        
                    elseif (self.indiced_rho_ini(block,3) == 1)
                        
                        for b = self.indiced_rho_ini(block,2):-1:self.indiced_rho_ini(block,1)
                            % For k>k_c case, From the last one
                            arraym = self.m_initial_con_cf(b,self.t_max_traj(k), self.x_max_traj(k));
                            array2 = self.subtractArray(arraym,temp_array);
                            if(~isempty(array2))
                                rows = size(list3,1);
                                list3(rows+1,:) = array2;
                                rows = size(tempM,1);
                                tempM(rows+1,:) = arraym;
                                break;
                            end
                        end
                        
                    elseif isnan(self.indiced_rho_ini(block,3))
                        
                        for b = self.indiced_rho_ini(block,1):self.indiced_rho_ini(block,2)
                            % For undefined initial blocks, check all
                            arraym = self.m_initial_con_ff(b,self.t_max_traj(k), self.x_max_traj(k));
                            array2 = self.subtractArray(arraym,temp_array);
                            if(~isempty(array2))
                                rows = size(list3,1);
                                list3(rows+1,:) = array2;
                                rows = size(tempM,1);
                                tempM(rows+1,:) = arraym;
                            end
                            
                            arraym = self.m_initial_con_cf(b,self.t_max_traj(k), self.x_max_traj(k));
                            array2 = self.subtractArray(arraym,temp_array);
                            if(~isempty(array2))
                                rows = size(list3,1);
                                list3(rows+1,:) = array2;
                                rows = size(tempM,1);
                                tempM(rows+1,:) = arraym;
                            end
                        end
                        
                    end
                    
                end
                

                % Solutions associated with downstream conditions
                for n=1:self.num_ds_con
                    arraym = self.m_ds_con(n,self.t_max_traj(k), self.x_max_traj(k));
                    array2 = self.subtractArray(arraym,temp_array);
                    if(~isempty(array2) && arraym(1, self.dv_link.downstream(1)-1 + n) ~= self.T_ds(n))
                        %Add a constraint
                        rows = size(list3,1);
                        list3(rows+1,:) = array2;
                        %Assign the solution to temporary Matrix
                        rows = size(tempM,1);
                        tempM(rows+1,:) = arraym;
                        break;
                    end
                end
                
                % Solutions associated with upstream conditions
                for n=1:self.num_us_con
                    arraym = self.m_us_con(n,self.t_max_traj(k), self.x_max_traj(k));
                    array2 = self.subtractArray(arraym, temp_array);
                    if(~isempty(array2)  && arraym(1,self.size_row)==0)
                        %Add a constraint
                        rows = size(list3,1);
                        list3(rows+1,:) = array2;
                        %Assign the solution to temporary Matrix
                        rows = size(tempM,1);
                        tempM(rows+1,:) = arraym;
                        break;
                    end
                    
                end
                
                % Solutions associated with internal trajectory conditions
                for n=1:self.num_internal_con
                    
                    if(n~=k)
                        arraym = self.m_traj_con(n,self.t_max_traj(k), self.x_max_traj(k));
                        array2 = self.subtractArray(arraym,temp_array);
                        if(~isempty(array2))
                            %Add a constraint
                            rows = size(list3,1);
                            list3(rows+1,:) = array2;
                            %Assign the solution to temporary Matrix
                            rows = size(tempM,1);
                            tempM(rows+1,:) = arraym;
                            
                        end
                        
                    end
                    
                    
                end
                
                % Solutions associated with density conditions
                index_u = find(self.t_dens ~= self.end_time);
                for i=1:length(index_u) % For density
                    
                    n = index_u(i);
                    
                    arraym = self.m_dens_con_ff(n,self.t_max_traj(k), self.x_max_traj(k));
                    array2 = self.subtractArray(arraym,temp_array);
                    if(~isempty(array2))
                        rows = size(list3,1);
                        list3(rows+1,:) = array2;
                        rows = size(tempM,1);
                        tempM(rows+1,:) = arraym;
                    end
                    
                    arraym = self.m_dens_con_cf(n,self.t_max_traj(k), self.x_max_traj(k));
                    array2 = self.subtractArray(arraym,temp_array);
                    if(~isempty(array2))
                        rows = size(list3,1);
                        list3(rows+1,:) = array2;
                        rows = size(tempM,1);
                        tempM(rows+1,:) = arraym;
                    end
                    
                end
                

                
                % Define all the possible combinations according to the solutions that apply
                % at the specified point
                % Lower bound constraints
                % Each row of omb_mat select one solution, for instance (b1, b2) = (1, 0)
                % (1-b1)*inf + b2*inf + L >= M2 selects M2 solution.         
                rowsM = size(tempM,1);
                for i = 1:rowsM
                    
                    array = temp_array;
                    
                    %Decode the binary combinations
                    for counter=1:self.nb_max_traj(1,k)
                        if comb_mat(i,nb-self.nb_max_traj(1,k)+counter) == 1
                            array(1, self.dv_link.bool_internal_min(1)-1 + countB + counter) = -Cmax;
                            
                        elseif comb_mat(i,nb-self.nb_max_traj(1,k)+counter) == 0
                            array(1, self.dv_link.bool_internal_min(1)-1 + countB + counter) = Cmax;
                            
                        end
                    end

                    array(1,self.size_row) = -Cmax*(sum(comb_mat(i,:))); %RHS (set negative to be on same side)
                    % Add constraint to the matrix constraints
                    array2 = self.subtractArray(array,tempM(i,:));
                    rows = size(list3,1);
                    list3(rows+1,:) = array2;
                end
                
                %Define the last constraint to bound the binary combination
                array = zeros(1,self.size_row);
                % Define constraint matrix elements for binary terms
                for counter=1:self.nb_max_traj(1,k)
                    
                    array(1, self.dv_link.bool_internal_min(1)-1 + countB + counter) = -2^(self.nb_max_traj(1,k)-counter);
                    
                end
                array(1,self.size_row) = -(rowsM-1); %RHS
                rows = size(list3,1);
                list3(rows+1,:) = array;
                
                countB = countB + self.nb_max_traj(k);
            end

        end
        


        %==========================================================================
        % Function to create the MILP constraints associated with density
        % conditions. Each density condition end point has a
        % vehicle id L.
        % L = min( Mi (all value condition solutions) )
        % Use binary variables to remove nonliear min function as
        %   L <= Mi for all i
        %   inf*bi + L >= Mi for all i
        %   sum(1-bi) = 1 (one and only one bi = 0)
        % output:
        %       list[:,1:end-1] * x >= list[:,end]
        function [list3] = setDensityConstraints(self)
            
            list3 = zeros(0,0);
            
            %Binary variables useful counter as a reference
            countB = 0;
            Cmax = inf; 
            
            % nb: number of binary variables per point
            nb = ceil(log2(2*(self.num_initial_con+1) + ...
                    (self.num_us_con+1)+ (self.num_ds_con+1) +...
                    (self.num_internal_con+1)+(self.num_density_con+1)));
            
            % comb_mat: combination matrix for possible binary combinations
            comb_mat = zeros(2^nb,nb);
            
            for i = 1:2^nb
                comb_mat(i,:) = double(dec2bin(i-1,nb)=='1');
            end
            
            % Binary inequalities induced by density conditions
            % For each density condition
            for k = 1:self.num_density_con
                
                 % ------------------------------------------------------------------------
                 %x_min_u and t_u point------------------
                if (k==1 || (self.x_min_dens(k)~=self.x_max_dens(k-1) &&...
                        self.t_dens(k)~=self.t_dens(k-1)))
                    tempM = zeros(0,0);
                    %Upper bound constraints
                    %L1<=M1
                    %L1<=M2
                    temp_array = self.dens_con(k, self.t_dens(k), self.x_min_dens(k));
                    
                    % Initial condition solutions Mi
                    % Lk <= Mi
                    for block = 1:size(self.indiced_rho_ini,1) 
                        
                        if (self.indiced_rho_ini(block,3) == 0)
                            
                            for b = self.indiced_rho_ini(block,1):self.indiced_rho_ini(block,2)
                                %-1 to be consistent with other functions self.mtau...
                                arraym = self.m_initial_con_ff(b,self.t_dens(k), self.x_min_dens(k));
                                array2 = self.subtractArray(arraym,temp_array);
                                if(~isempty(array2))
                                    rows = size(list3,1);
                                    list3(rows+1,:) = array2;
                                    rows = size(tempM,1);
                                    tempM(rows+1,:) = arraym;
                                    break;  % Only need to consider the first one
                                end
                            end
                            
                        elseif (self.indiced_rho_ini(block,3) == 1)
                            
                            for b = self.indiced_rho_ini(block,2):-1:self.indiced_rho_ini(block,1)
                                % For k>k_c case, From the last one
                                arraym = self.m_initial_con_cf(b,self.t_dens(k), self.x_min_dens(k));
                                array2 = self.subtractArray(arraym,temp_array);
                                if(~isempty(array2))
                                    rows = size(list3,1);
                                    list3(rows+1,:) = array2;
                                    rows = size(tempM,1);
                                    tempM(rows+1,:) = arraym;
                                    break;
                                end
                            end
                            
                        elseif (isnan(self.indiced_rho_ini(block,3)))
                            
                            for b = self.indiced_rho_ini(block,1):self.indiced_rho_ini(block,2)
                                % For undefined initial blocks, check all
                                arraym = self.m_initial_con_ff(b,self.t_dens(k), self.x_min_dens(k));
                                array2 = self.subtractArray(arraym,temp_array);
                                if(~isempty(array2))
                                    rows = size(list3,1);
                                    list3(rows+1,:) = array2;
                                    rows = size(tempM,1);
                                    tempM(rows+1,:) = arraym;
                                end
                                
                                arraym = self.m_initial_con_cf(b,self.t_dens(k), self.x_min_dens(k));
                                array2 = self.subtractArray(arraym,temp_array);
                                if(~isempty(array2))
                                    rows = size(list3,1);
                                    list3(rows+1,:) = array2;
                                    rows = size(tempM,1);
                                    tempM(rows+1,:) = arraym;
                                end
                            end
                            
                        end
                    end 
                    

                    % Downstream condition solutions Mi
                    % Lk <= Mi
                    for n=1:self.num_ds_con    
                        
                        arraym = self.m_ds_con(n,self.t_dens(k), self.x_min_dens(k));
                        array2 = self.subtractArray(arraym,temp_array);
                        if(~isempty(array2) && arraym(1,self.dv_link.downstream(1)-1 + n) ~= self.T_ds(n))
                            %Add a constraint
                            rows = size(list3,1);
                            list3(rows+1,:) = array2;
                            %Assign the solution to temporary Matrix
                            rows = size(tempM,1);
                            tempM(rows+1,:) = arraym;
                            break;  %only one
                        end
                        
                    end
                    

                    % Upstream condition solutions Mi
                    % Lk <= Mi
                    for n=1:self.num_us_con
                        
                        arraym = self.m_us_con(n,self.t_dens(k), self.x_min_dens(k));
                        array2 = self.subtractArray(arraym, temp_array);
                        if(~isempty(array2) && arraym(1,self.size_row)==0)
                            %Add a constraint
                            rows = size(list3,1);
                            list3(rows+1,:) = array2;
                            %Assign the solution to temporary Matrix
                            rows = size(tempM,1);
                            tempM(rows+1,:) = arraym;
                            break;
                        end
                        
                    end
                    

                    % Internal condition solutions Mi
                    % Lk <= Mi
                    for n=1:self.num_internal_con    
                        
                        arraym = self.m_traj_con(n,self.t_dens(k), self.x_min_dens(k));
                        array2 = self.subtractArray(arraym,temp_array);
                        if(~isempty(array2))
                            rows = size(list3,1);
                            list3(rows+1,:) = array2;
                            rows = size(tempM,1);
                            tempM(rows+1,:) = arraym;
                        end
                        
                    end
                    

                    % Density condition solutions Mi
                    % Lk <= Mi
                    % Find density conditions that are not at t_end. Not need to compute solutions
                    % associated with density conditions at t_end
                    index_u = find(self.t_dens <= self.end_time);
                    for i=1:length(index_u) % For density
                
                        n = index_u(i);
                        
                        
                        arraym = self.m_dens_con_ff(n,self.t_dens(k), self.x_min_dens(k));
                        array2 = self.subtractArray(arraym,temp_array);
                        if(~isempty(array2))
                            rows = size(list3,1);
                            list3(rows+1,:) = array2;
                            rows = size(tempM,1);
                            tempM(rows+1,:) = arraym;
                        end
                            
                        arraym = self.m_dens_con_cf(n,self.t_dens(k), self.x_min_dens(k));
                        array2 = self.subtractArray(arraym,temp_array);
                        if(~isempty(array2))
                            rows = size(list3,1);
                            list3(rows+1,:) = array2;
                            rows = size(tempM,1);
                            tempM(rows+1,:) = arraym;
                        end
                       
                        
                    end
                    
                    
                    
                    % Define all the possible combinations according to the solutions that apply
                    % at the specified point
                    % Lower bound constraints
                    % Each row of omb_mat select one solution, for instance (b1, b2) = (1, 0)
                    % (1-b1)*inf + b2*inf + L >= M2 selects M2 solution.       
                    rowsM = size(tempM,1);
                    for i = 1:rowsM
                        
                        array = temp_array;
                        
                        %Decode the binary combinations
                        for counter=1:self.nb_min_dens(k)
                            if comb_mat(i,nb-self.nb_min_dens(k+1)+counter) == 1
                                array(1, self.dv_link.bool_density_min(1)-1 + countB + counter) = -Cmax;
                                
                            elseif comb_mat(i,nb-self.nb_min_dens(k+1)+counter) == 0
                                array(1, self.dv_link.bool_density_min(1)-1 + countB+ counter) = Cmax;
                                
                            end
                        end

                        array(1,self.size_row) = -Cmax*(sum(comb_mat(i,:))); %RHS (set negative to be on same side)
                        % Add constraint to the MILP matrix
                        array2 = self.subtractArray(array,tempM(i,:));
                        rows = size(list3,1);
                        list3(rows+1,:) = array2;
                    end
                    
                    %Define the last constraint to bound the binary combination
                    array = zeros(1,self.size_row);
                    
                    for counter=1:self.nb_min_dens(k)
                        
                        array(1, self.dv_link.bool_density_min(1)-1 + countB + counter) = ...
                            -2^(self.nb_min_dens(k)-counter);
                        
                    end
                    
                    array(1,self.size_row) = -(rowsM-1); %RHS (maximum possible value of the binary comb)
                    rows = size(list3,1);
                    list3(rows+1,:) = array;
                    
                    countB = countB+self.nb_min_dens(k);
                    
                end


                % ------------------------------------------------------------------------
                %x_max_u and t_u point-------------------------
                tempM = zeros(0,0);
                
                %Upper bound constraints
                %L1<=M1
                %L1<=M2
                temp_array = self.dens_con(k, self.t_dens(k), self.x_max_dens(k));
                

                % Solutions associated with initial conditions
                for block = 1:size(self.indiced_rho_ini,1) 
                    
                    if (self.indiced_rho_ini(block,3) == 0)
                        
                        for b = self.indiced_rho_ini(block,1):self.indiced_rho_ini(block,2)
                            %-1 to be consistent with other functions self.mtau...
                            arraym = self.m_initial_con_ff(b,self.t_dens(k), self.x_max_dens(k));
                            array2 = self.subtractArray(arraym,temp_array);
                            if(~isempty(array2))
                                rows = size(list3,1);
                                list3(rows+1,:) = array2;
                                rows = size(tempM,1);
                                tempM(rows+1,:) = arraym;
                                break;  % Only need to consider the first one
                            end
                        end
                        
                    elseif (self.indiced_rho_ini(block,3) == 1)
                        
                        for b = self.indiced_rho_ini(block,2):-1:self.indiced_rho_ini(block,1)
                            % For k>k_c case, From the last one
                            arraym = self.m_initial_con_cf(b,self.t_dens(k), self.x_max_dens(k));
                            array2 = self.subtractArray(arraym,temp_array);
                            if(~isempty(array2))
                                rows = size(list3,1);
                                list3(rows+1,:) = array2;
                                rows = size(tempM,1);
                                tempM(rows+1,:) = arraym;
                                break;
                            end
                        end
                        
                    elseif isnan(self.indiced_rho_ini(block,3))
                        
                        for b = self.indiced_rho_ini(block,1):self.indiced_rho_ini(block,2)
                            % For undefined initial blocks, check all
                            arraym = self.m_initial_con_ff(b,self.t_dens(k), self.x_max_dens(k));
                            array2 = self.subtractArray(arraym,temp_array);
                            if(~isempty(array2))
                                rows = size(list3,1);
                                list3(rows+1,:) = array2;
                                rows = size(tempM,1);
                                tempM(rows+1,:) = arraym;
                            end
                            
                            arraym = self.m_initial_con_cf(b,self.t_dens(k), self.x_max_dens(k));
                            array2 = self.subtractArray(arraym,temp_array);
                            if(~isempty(array2))
                                rows = size(list3,1);
                                list3(rows+1,:) = array2;
                                rows = size(tempM,1);
                                tempM(rows+1,:) = arraym;
                            end
                        end
                        
                    end
                    
                end
                
                % Solutions associated with downstream conditions
                for n=1:self.num_ds_con
                    arraym = self.m_ds_con(n,self.t_dens(k), self.x_max_dens(k));
                    array2 = self.subtractArray(arraym,temp_array);
                    if(~isempty(array2) && array(1,self.dv_link.downstream(1)-1 + n) ~= self.T_ds(n))
                        %Add a constraint
                        rows = size(list3,1);
                        list3(rows+1,:) = array2;
                        %Assign the solution to temporary Matrix
                        rows = size(tempM,1);
                        tempM(rows+1,:) = arraym;
                        break;
                    end
                end
                

                % Solutions associated with upstream conditions
                for n=1:self.num_us_con
                    arraym = self.m_us_con(n,self.t_dens(k), self.x_max_dens(k));
                    array2 = self.subtractArray(arraym, temp_array);
                    if(~isempty(array2)  && arraym(1,self.size_row)==0)
                        %Add a constraint
                        rows = size(list3,1);
                        list3(rows+1,:) = array2;
                        %Assign the solution to temporary Matrix
                        rows = size(tempM,1);
                        tempM(rows+1,:) = arraym;
                        break;
                    end
                    
                end
                
                 % Solutions associated with internal condition
                for n=1:self.num_internal_con
                    
                    arraym = self.m_traj_con(n,self.t_dens(k), self.x_max_dens(k));
                    array2 = self.subtractArray(arraym,temp_array);
                    if(~isempty(array2))
                        %Add a constraint
                        rows = size(list3,1);
                        list3(rows+1,:) = array2;
                        %Assign the solution to temporary Matrix
                        rows = size(tempM,1);
                        tempM(rows+1,:) = arraym;
                    end
                    
                end
                

                 % Solutions associated with density conditions
                index_u = find(self.t_dens <= self.end_time);
                for i=1:length(index_u) % For density
                    
                    n = index_u(i);
                    
                        arraym = self.m_dens_con_ff(n,self.t_dens(k), self.x_max_dens(k));
                        array2 = self.subtractArray(arraym,temp_array);
                        if(~isempty(array2))
                            rows = size(list3,1);
                            list3(rows+1,:) = array2;
                            rows = size(tempM,1);
                            tempM(rows+1,:) = arraym;
                        end
                        
                        arraym = self.m_dens_con_cf(n,self.t_dens(k), self.x_max_dens(k));
                        array2 = self.subtractArray(arraym,temp_array);
                        if(~isempty(array2))
                            rows = size(list3,1);
                            list3(rows+1,:) = array2;
                            rows = size(tempM,1);
                            tempM(rows+1,:) = arraym;
                        end
                 
                    
                end
                

                % Define all the possible combinations according to the solutions that apply
                % at the specified point
                % Lower bound constraints
                % Each row of omb_mat select one solution, for instance (b1, b2) = (1, 0)
                % (1-b1)*inf + b2*inf + L >= M2 selects M2 solution.     
                rowsM = size(tempM,1);
                for i = 1:rowsM
                    
                    array = temp_array;
                    
                    %Decode the binary combinations
                    % still use self.dv_link.bool_density_min(1) is because
                    % countB is not reset
                    for counter=1:self.nb_max_dens(1,k)
                        if comb_mat(i,nb-self.nb_max_dens(1,k)+counter) == 1
                            array(1, self.dv_link.bool_density_min(1)-1 + countB + counter) = -Cmax;
                            
                        elseif comb_mat(i,nb-self.nb_max_dens(1,k)+counter) == 0
                            array(1, self.dv_link.bool_density_min(1)-1 + countB + counter) = Cmax;
                            
                        end
                    end
                    array(1,self.size_row) = -Cmax*(sum(comb_mat(i,:))); %RHS (set negative to be on same side)
                    
                    % Add constraint to the matrix constraints
                    
                    array2 = self.subtractArray(array,tempM(i,:));
                    rows = size(list3,1);
                    list3(rows+1,:) = array2;
                end
                
                %Define the last constraint to bound the binary combination
                
                array = zeros(1,self.size_row);
                
                % Define constraint matrix elements for binary terms
                for counter=1:self.nb_max_dens(1,k)
                    
                    array(1, self.dv_link.bool_density_min(1)-1 + countB + counter) = ...
                         -2^(self.nb_max_dens(1,k)-counter);
                    
                end
                
                array(1,self.size_row) = -(rowsM-1); %RHS
                rows = size(list3,1);
                list3(rows+1,:) = array;
                
                countB = countB + self.nb_max_dens(k);
          
            end
            
        end
        
        
        
        %==========================================================================
        % Function to create hard constraitns for queue limit.
        % Intuition: if the queue extends to the point, then the solution L
        % at that point comes from the downstream. Since by triangular FD,
        % the point will at only in the charicateristic domain M_us of one value
        % condition. Similarly, it locates only in the charactheristic
        % domain M_ds of one downstream condition. Then we only need
        %       M_ds >= M_us
        % NOTE: this only works with only initial and boundary conditions
        % input: 
        %       max_length: the maximal length of the queue in meters
        %       from the downstream
        % output:
        %       the constraint matrix
        %       list[:,1:end-1] * x >= list[:,end]
        function [list] = setHardQueueLimit(self, max_length)
            
            len_link = self.ds_pos_m - self.us_pos_m;
            
            % create a set of queue limit point evenly discretized in 30 s 
            % in time
            T_queue_cum = self.start_time:30:self.end_time;
            if T_queue_cum(end) ~= self.end_time
                T_queue_cum = [T_queue_cum'; self.end_time];
            end
            num_queue_pt = length(T_queue_cum)-1;
            t_queue = T_queue_cum(2:end);
            
            if max_length >= len_link
                warning('ineffective queue limit')
            else
                x_queue = self.ds_pos_m - ones(num_queue_pt,1)*max_length;
            end
            
            list = zeros(num_queue_pt, self.size_row);
            rows = 0;
            
            % for each point
            for pt = 1: num_queue_pt
                
                % Note only one of M_ds and M_init_cong will be non empty
                % and only one of M_us and M_init_free will be non empty
                
                % Downstream characteristic domain
                M_ds = [];
                for ds = 1:self.num_ds_con
                    tmp_ds = self.m_ds_con(ds, t_queue(pt), x_queue(pt));
                    if ~isempty(tmp_ds) && tmp_ds(1, self.dv_link.downstream(1)-1 + ds) ~= self.T_ds(ds)
                        M_ds = tmp_ds;
                        break;  % Got the solution in the characteristic domain
                    end
                end
                
                % If not in downstream characteristic domain. Then must be
                % in the initial condtion domain when its congested. 
                if isempty(M_ds)
                    % congested initial conditions;
                    M_init_cong = [];
                    for b = 1:self.num_initial_con
                        tmp_init_cong = self.m_initial_con_cf(b, t_queue(pt), x_queue(pt));
                        if ~isempty(tmp_init_cong) && ...
                                tmp_init_cong(1, self.dv_link.initial(1)-1 + b) ~= -self.X_grid(b)
                            M_init_cong = tmp_init_cong;
                            break;
                        end
                    end
                    
                    % impossible. something must be wrong
                    if isempty(M_init_cong)
                        error('Check the queue limit constraints.')
                    else
                        M_ds = M_init_cong;
                    end
                    
                end
                
                % Upstream characteristic domain
                M_us = [];
                for us = 1:self.num_us_con
                    M_us = self.m_us_con(us, t_queue(pt), x_queue(pt));
                    if ~isempty(M_us) && M_us(1, self.size_row) == 0
                        break;  % Got the solution in the characteristic domain
                    end
                end
                
                if isempty(M_us)
                    % freeflow initial conditions
                    M_init_free = [];
                    for b = 1:self.num_initial_con
                        tmp_init_free = self.m_initial_con_ff(b, t_queue(pt), x_queue(pt));
                        if ~isempty(tmp_init_free) && ...
                                tmp_init_free(1, self.size_row) == 0
                            M_init_free = tmp_init_free;
                            break;
                        end
                    end
                    
                    % impossible. something must be wrong
                    if isempty(M_init_free)
                        error('Check the queue limit constraints.')
                    else
                        M_us = M_init_free;
                    end
                    
                end
                
                % Note the sign: M_ds >= M_us
                % list(:,1:end-1) >= list(:, end)
                array = self.subtractArray(M_ds, M_us);
                
                rows = rows+1;
                list(rows,:) = array; 
                
            end
            
            % truncate matrix
            list = list(1:rows,:);
            
            
        end


        
        %==========================================================================
        % Function to create soft constraints for queue limits
        % For any point (x,t), define vehicle label L, and slack variable s
        % constraints: L <= M^{us(+initial_free), ds(+initial_cong) }; 
        %              L = M^us - s; s>=0; 
        % objective: min s
        % input:
        % output: 
        %       list: the constraint matrix; 
        %       list[:,1:end-1] * x >= list[:,end]
        %       weight: the weight for q1, this is used for constructing
        %           the entropy condition
        function [list, weight] = setSoftQueueLimit(self)
            
            % Each point introduces 5 inequalities
            list = zeros(self.num_soft_queue_pt*5, self.size_row);
            rows = 0;
            % keep track of the weight of q1 which is in M^us
            list_us = zeros(self.num_soft_queue_pt, self.size_row);
            rows_us = 0;
                        
            % for each point
            for pt = 1:self.num_soft_queue_pt
                
                % only consider the points in the future.
                if self.soft_queue_t(pt) > self.now_time
                    
                    % Note only one of M_ds and M_init_cong will be non empty
                    % and only one of M_us and M_init_free will be non empty
                    
                    % Downstream characteristic domain
                    M_ds = [];
                    for ds = 1:self.num_ds_con
                        tmp_ds = self.m_ds_con(ds, self.soft_queue_t(pt), self.soft_queue_x(pt));
                        if ~isempty(tmp_ds) && ...
                                tmp_ds(1, self.dv_link.downstream(1)-1 + ds) ~= self.T_ds(ds)
                            M_ds = tmp_ds;
                            break;  % Got the solution in the characteristic domain
                        end
                    end
                    
                    % If not in downstream characteristic domain. Then must be
                    % in the initial condtion domain when its congested.
                    if isempty(M_ds)
                        % congested initial conditions;
                        M_init_cong = [];
                        for b = 1:self.num_initial_con
                            tmp_init_cong = self.m_initial_con_cf(b, self.soft_queue_t(pt), self.soft_queue_x(pt));
                            if ~isempty(tmp_init_cong) && ...
                                    tmp_init_cong(1, self.dv_link.initial(1)-1 + b) ~= -self.X_grid(b)
                                M_init_cong = tmp_init_cong;
                                break;
                            end
                        end
                        
                        % impossible. something must be wrong
                        if isempty(M_init_cong)
                            error('Check the queue limit constraints.')
                        else
                            M_ds = M_init_cong;
                        end
                        
                    end
                    
                    % Upstream characteristic domain
                    M_us = [];
                    for us = 1:self.num_us_con
                        M_us = self.m_us_con(us, self.soft_queue_t(pt), self.soft_queue_x(pt));
                        if ~isempty(M_us) && M_us(1, self.size_row) == 0
                            break;  % Got the solution in the characteristic domain
                        end
                    end
                    
                    if isempty(M_us)
                        % freeflow initial conditions
                        M_init_free = [];
                        for b = 1:self.num_initial_con
                            tmp_init_free = self.m_initial_con_ff(b, self.soft_queue_t(pt), self.soft_queue_x(pt));
                            if ~isempty(tmp_init_free) && ...
                                    tmp_init_free(1, self.size_row) == 0
                                M_init_free = tmp_init_free;
                                break;
                            end
                        end
                        
                        % impossible. something must be wrong
                        if isempty(M_init_free)
                            error('Check the queue limit constraints.')
                        else
                            M_us = M_init_free;
                        end
                        
                    end
                    
                    %=========add constraints==========
                    % M_us - L >= 0
                    array = M_us;
                    array(1, self.dv_link.queue_L(1)-1 + pt) = -1;
                    rows = rows+1;
                    list(rows,:) = array;
                    % keep track of weights
                    rows_us = rows_us+1;
                    list_us(rows_us,:) = array;
                    
                    % M_ds - L >= 0
                    array = M_ds;
                    array(1, self.dv_link.queue_L(1)-1 + pt) = -1;
                    rows = rows+1;
                    list(rows,:) = array;
                    
                    % M_us - s - L >= & L - M_us + s >=0
                    array = M_us;
                    array(1, self.dv_link.queue_L(1)-1 + pt) = -1;
                    array(1, self.dv_link.queue_s(1)-1 + pt) = -1;
                    rows = rows+1;
                    list(rows,:) = array;
                    array = -array;
                    rows = rows+1;
                    list(rows,:) = array;
                    
                    % s >= 0
                    array = zeros(1, self.size_row);
                    array(1, self.dv_link.queue_s(1)-1 + pt) = 1;
                    rows = rows+1;
                    list(rows,:) = array;
                    
                end
                
            end
            
            % truncate matrix
            list = list(1:rows,:);
            list_us = list_us(1:rows_us,:);
            
            % compute the sum of the weight
            if ~isempty(list_us)
                weight = sum( list_us(:,...
                    self.dv_link.upstream(1):...
                    self.dv_link.upstream(2)) ,1);
                weight = weight';
            else
                weight = [];
            end
            
        end
    
        
        
        
        %=========================================================================
        %Get Number Binary variables
        
        function [nb_min,nb_max, nb_min_u, nb_max_u] = getBinaryvar(self)
            
            nb_min = zeros(1,self.num_internal_con);
            nb_max = zeros(1,self.num_internal_con);
            
            nb_min_u = zeros(1,self.num_density_con);
            nb_max_u = zeros(1,self.num_density_con);
            
            % Binaries induced by internal conditions
            for k = 1:self.num_internal_con
                % For x_min
                if (k==1 || (self.x_min_traj(k)~=self.x_max_traj(k-1) &&...
                        self.t_min_traj(k)~=self.t_max_traj(k-1)))
                    
                    countM=0;
                    
                    % For initial
                    for block = 1:size(self.indiced_rho_ini,1) 
                        
                        if (self.indiced_rho_ini(block,3) == 0)   % For k<k_c case, only the first one
                            
                            for b = self.indiced_rho_ini(block,1):self.indiced_rho_ini(block,2)
                                array = self.m_initial_con_ff(b,self.t_min_traj(k), self.x_min_traj(k));
                                if(~isempty(array))
                                    countM = countM+1;
                                    break;  
                                end
                            end
                            
                        elseif (self.indiced_rho_ini(block,3) == 1) % For k>k_c case, only the last one
                            
                            for b = self.indiced_rho_ini(block,2):-1:self.indiced_rho_ini(block,1)
                                array = self.m_initial_con_cf(b,self.t_min_traj(k), self.x_min_traj(k));
                                if(~isempty(array))
                                    countM = countM+1;
                                    break;
                                end
                            end
                            
                        elseif isnan(self.indiced_rho_ini(block,3))
                            
                            for b = self.indiced_rho_ini(block,1):self.indiced_rho_ini(block,2)
                                array = self.m_initial_con_ff(b,self.t_min_traj(k), self.x_min_traj(k));
                                if(~isempty(array))
                                    countM = countM+1;
                                end
                                
                                array = self.m_initial_con_cf(b,self.t_min_traj(k), self.x_min_traj(k));
                                if(~isempty(array))
                                    countM = countM+1;
                                end
                            end
                            
                        end
                        
                    end
                    
                    % for downstream
                    % in downstream characteristic domain
                    for n=1:self.num_ds_con
                        
                        array = self.m_ds_con(n,self.t_min_traj(k), self.x_min_traj(k));
                        if(~isempty(array) && array(1, self.dv_link.downstream(1)-1 + n) ~= self.T_ds(n))
                            % If in characteristic domain
                            countM = countM+1;
                            break;
                        end
                        
                    end
                    
                    % for upstream
                    % in upstream characteristic domain
                    for n=1:self.num_us_con
                        
                        array = self.m_us_con(n,self.t_min_traj(k), self.x_min_traj(k));
                        if(~isempty(array) && array(1,n) ~= self.T_us(n) )
                            countM = countM+1;
                            break;
                        end
                        
                    end
                    
                    % for other internal conditions
                    for n=1:self.num_internal_con
                        
                        if (n~=k)
                            array = self.m_traj_con(n,self.t_min_traj(k), self.x_min_traj(k));
                            if(~isempty(array))
                                countM = countM+1;
                            end
                        end
                        
                    end
                    
                    % for density conditions
                    index_u = find(self.t_dens~= self.end_time);
                    for i=1:length(index_u) % For density
                
                        n = index_u(i);
                        
                        array = self.m_dens_con_ff(n,self.t_min_traj(k), self.x_min_traj(k));
                        if(~isempty(array))
                            countM = countM+1;
                        end
                        
                        array = self.m_dens_con_cf(n,self.t_min_traj(k), self.x_min_traj(k));
                        if(~isempty(array))
                            countM = countM+1;
                        end
                        
                    end
                    
                    nb_min(1,k) = ceil(log2(countM));
                    
                end
                
                countM=0;
                
                %Upper bound constraints
                %L1<=M1
                %L1<=M2
                
                % For initial
                for block = 1:size(self.indiced_rho_ini,1)
                    
                    if (self.indiced_rho_ini(block,3) == 0)   % For k<k_c case, only the first one
                        
                        for b = self.indiced_rho_ini(block,1):self.indiced_rho_ini(block,2)
                            array = self.m_initial_con_ff(b,self.t_min_traj(k), self.x_min_traj(k));
                            if(~isempty(array))
                                countM = countM+1;
                                break;
                            end
                        end
                        
                    elseif (self.indiced_rho_ini(block,3) == 1) % For k>k_c case, only the last one
                        
                        for b = self.indiced_rho_ini(block,2):-1:self.indiced_rho_ini(block,1)
                            array = self.m_initial_con_cf(b,self.t_min_traj(k), self.x_min_traj(k));
                            if(~isempty(array))
                                countM = countM+1;
                                break;
                            end
                        end
                        
                    elseif isnan(self.indiced_rho_ini(block,3))
                        
                        for b = self.indiced_rho_ini(block,1):self.indiced_rho_ini(block,2)
                            array = self.m_initial_con_ff(b,self.t_min_traj(k), self.x_min_traj(k));
                            if(~isempty(array))
                                countM = countM+1;
                            end
                            
                            array = self.m_initial_con_cf(b,self.t_min_traj(k), self.x_min_traj(k));
                            if(~isempty(array))
                                countM = countM+1;
                            end
                        end
                        
                    end
                    
                end
                
                % for downstream
                % in downstream characteristic domain
                for n=1:self.num_ds_con
                    
                    array = self.m_ds_con(n,self.t_min_traj(k), self.x_min_traj(k));
                    if(~isempty(array) && array(1,self.dv_link.downstream(1)-1 + n) ~= self.T_ds(n))
                        % If in characteristic domain
                        countM = countM+1;
                        break;
                    end
                    
                end
                
                % for upstream
                % in upstream characteristic domain
                for n=1:self.num_us_con
                    
                    array = self.m_us_con(n,self.t_min_traj(k), self.x_min_traj(k));
                    if(~isempty(array) && array(1,n) ~= self.T_us(n) )
                        countM = countM+1;
                        break;
                    end
                    
                end
                
                % for other internal conditions
                for n=0:self.num_internal_con
                    
                    if (n~=k)
                        array = self.m_traj_con(n,self.t_min_traj(k), self.x_min_traj(k));
                        if(~isempty(array))
                            countM = countM+1;
                        end
                    end
                    
                end
                
                % for density conditions
                index_u = find(self.t_dens~= self.end_time);
                for i=1:length(index_u) % For density
                    
                    n = index_u(i);
                    
                    array = self.m_dens_con_ff(n,self.t_min_traj(k), self.x_min_traj(k));
                    if(~isempty(array))
                        countM = countM+1;
                    end
                    
                    array = self.m_dens_con_cf(n,self.t_min_traj(k), self.x_min_traj(k));
                    if(~isempty(array))
                        countM = countM+1;
                    end
                    
                end
                
                nb_max(1,k) = ceil(log2(countM));
            end
            
            
            %Binaries induced by density conditions
            for k = 1:self.num_density_con
                % For x_min
                if (k==1 || (self.x_min_dens(k)~=self.x_max_dens(k-1) && self.t_dens(k)~=self.t_dens(k-1)))
                    
                    countM=0;
                    
                    for block = 1:size(self.indiced_rho_ini,1) % For initial
                        
                        if (self.indiced_rho_ini(block,3) == 0)   % For k<k_c case, only the first one
                            
                            for b = self.indiced_rho_ini(block,1):self.indiced_rho_ini(block,2)
                                array = self.m_initial_con_ff(b,self.t_dens(k), self.x_min_dens(k));
                                if(~isempty(array))
                                    countM = countM+1;
                                    break;  
                                end
                            end
                            
                        elseif (self.indiced_rho_ini(block,3) == 1) % For k>k_c case, only the last one
                            
                            for b = self.indiced_rho_ini(block,2):-1:self.indiced_rho_ini(block,1)
                                array = self.m_initial_con_cf(b,self.t_dens(k), self.x_min_dens(k));
                                if(~isempty(array))
                                    countM = countM+1;
                                    break;
                                end
                            end
                            
                        elseif isnan(self.indiced_rho_ini(block,3))
                            
                            for b = self.indiced_rho_ini(block,1):self.indiced_rho_ini(block,2)
                                array = self.m_initial_con_ff(b,self.t_dens(k), self.x_min_dens(k));
                                if(~isempty(array))
                                    countM = countM+1;
                                end
                                
                                array = self.m_initial_con_cf(b,self.t_dens(k), self.x_min_dens(k));
                                if(~isempty(array))
                                    countM = countM+1;
                                end
                            end
                            
                        end
                        
                    end
                    
                    % downstream
                    for n=1:self.num_ds_con
                        
                        array = self.m_ds_con(n,self.t_dens(k), self.x_min_dens(k));
                        if(~isempty(array) && array(1,self.dv_link.downstream(1)-1 + n) ~= self.T_ds(n))
                            countM = countM+1;
                            break;
                        end
                        
                    end
                    
                    % upstream
                    for n=1:self.num_us_con
                        
                        array = self.m_us_con(n,self.t_dens(k), self.x_min_dens(k));
                        if(~isempty(array) && array(1,n) ~= self.T_us(n) )
                            countM = countM+1;
                            break;
                        end
                        
                    end
                    
                    % internal conditions
                    for n=1:self.num_internal_con
                        
                        array = self.m_traj_con(n,self.t_dens(k), self.x_min_dens(k));
                        if(~isempty(array))
                            countM = countM+1;
                        end
                        
                    end
                    
                    index_u = find(self.t_dens~= self.end_time);
                    for i=1:length(index_u) % For density
                
                        n = index_u(i);
                        
                  
                            array = self.m_dens_con_ff(n,self.t_dens(k), self.x_min_dens(k));
                            if(~isempty(array))
                                countM = countM+1;
                            end
                            
                            array = self.m_dens_con_cf(n,self.t_dens(k), self.x_min_dens(k));
                            if(~isempty(array))
                                countM = countM+1;
                            end
                  
                        
                    end
                    
                    nb_min_u(1,k) = ceil(log2(countM));
                    
                end
                
                countM=0;
                
                %Upper bound constraints
                %L1<=M1
                %L1<=M2
                
                for block = 1:size(self.indiced_rho_ini,1) % For initial
                    
                    if (self.indiced_rho_ini(block,3) == 0)
                        
                        for b = self.indiced_rho_ini(block,1):self.indiced_rho_ini(block,2)
                            array = self.m_initial_con_ff(b,self.t_dens(k), self.x_max_dens(k));
                            if(~isempty(array))
                                countM = countM+1;
                                break;  % Only need to consider the first one
                            end
                        end
                        
                    elseif (self.indiced_rho_ini(block,3) == 1)
                        
                        for b = self.indiced_rho_ini(block,2):-1:self.indiced_rho_ini(block,1)
                            % For k>k_c case, From the last one
                            array = self.m_initial_con_cf(b,self.t_dens(k), self.x_max_dens(k));
                            if(~isempty(array))
                                countM = countM+1;
                                break;
                            end
                        end
                        
                    elseif isnan(self.indiced_rho_ini(block,3))
                        
                        for b = self.indiced_rho_ini(block,1)-1:self.indiced_rho_ini(block,2)-1
                            array = self.m_initial_con_ff(b,self.t_dens(k), self.x_max_dens(k));
                            if(~isempty(array))
                                countM = countM+1;
                            end
                            
                            array = self.m_initial_con_cf(b,self.t_dens(k), self.x_max_dens(k));
                            if(~isempty(array))
                                countM = countM+1;
                            end
                        end
                        
                    end
                end
                
                % downstream
                for n=1:self.num_ds_con
                    array = self.m_ds_con(n,self.t_dens(k), self.x_max_dens(k));
                    if(~isempty(array) && array(1,self.dv_link.downstream(1)-1 + n) ~= self.T_ds(n))
                        countM = countM+1;
                        break;
                    end
                end
                
                % upstream
                for n=1:self.num_us_con
                    array = self.m_us_con(n,self.t_dens(k), self.x_max_dens(k));
                    if(~isempty(array) && array(1,n) ~= self.T_us(n) )
                        countM = countM+1;
                        break;
                    end
                end
                
                % internal
                for n=1:self.num_internal_con
                    
                        array = self.m_traj_con(n,self.t_dens(k), self.x_max_dens(k));
                        if(~isempty(array))
                            countM = countM+1;
                        end
              
                end
                
                index_u = find(self.t_dens ~= self.end_time);
                for i=1:length(index_u) % For density
                    
                    n = index_u(i);
                    
                 
                        array = self.m_dens_con_ff(n,self.t_dens(k), self.x_max_dens(k));
                        if(~isempty(array))
                            countM = countM+1;
                        end
                        
                        array = self.m_dens_con_cf(n,self.t_dens(k), self.x_max_dens(k));
                        if(~isempty(array))
                            countM = countM+1;
                        end
               
                    
                end
                
                nb_max_u(1,k) = ceil(log2(countM));
            end
      
        end
        

        
        %==========================================================================
        % Method to substract to arrays with a "null" condition
        function [output] = subtractArray(~,value1,value2)
            if ( (~isempty(value1)) && (~isempty(value2)) )
                output = value1-value2;
            else
                output = []; %return "null"
            end
        end
        
        
        
        %===============================================================
        % Utility function
        % This function is used to group the same element into blocks
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
        function indicedGroupValue = groupSameElement(~, k_ori,minlength)
            
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
                [values, ivalues, ~] = unique(k_tmp,'stable');
                
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
    
        
        
        
    end
end





