% Yanning Li, Sep 11, 2015

% This class sets the inequality constraints for the optimization program.
% Cumunlativley, it has the following features
% 1. It is based on the setIneqConstraints_grid.m which only adds the 
%   minimum set of constarints based on heuristic knowledge. See the code 
%   description in setIneqConstraints_grid.m for details.
% 2. The variable names are changed which makes the code much more
%   readable. The original code varialble names are defined in agreement
%   with the equations in paper. 
% 3. This version has full support for the initial condition, the boundary
%   condition, the internal condition, the density condition. The
%   horizontal traffic light condition was implemented in another version
%   of the code. 
% 4. This version supports different levels of error for the value
%   conditions.

% Sep 01, 2015
% Remark: this version does not support internal or density condition with
% uneven discretization of the boundary flows. It also assume all sensors
% have same level of accuracy.


classdef setIneqConstraints
    
    properties
        
        n_max_us;  % num of upstream steps 
        n_max_ds;  % num of downstream steps
        m_max;  % num of internal conditions
        b_max;  % num of initial conditions
        u_max;  % num of density conditions
        
        v;
        w;
        k_c;
        k_m;
        
        start_time;
        end_time;   % end of simulation time window
        
        xi;     % post meter of the upstream position
        chi;    % post meter of the downstream position
        
        % up and downstream time grid
        T_us;       %upstream_steps x 1
        T_us_cum;   % cumulative; upstream_steps+1 x 1
        T_ds;
        T_ds_cum;
        
        e_max;  % error for historical data, assume to be larger +-20%
        list;
        
        ini_segments;   % uneven space descritization of link b_max+1 x 1
        rho_ini;    % b_max x 1
        indiced_rho_ini;    
        X;
        
        % for boundary conditions
        qin_meas;
        qout_meas;
        
        % for internal conditions
        x_min;
        x_max;
        t_min;
        t_max;
        v_meas;
        r_meas;
        
        % for density conditions
        x_min_u;
        x_max_u;
        t_u;
        dens_meas;
        
        % auxiliary binary variables
        nb_min;
        nb_max;
        nb_min_u;
        nb_max_u;
        
        % auxiliary variables
        num_aux;
        size_row;
    end
    
    methods
        
        % This function set the inequality constraints of one link
        % input: 
        %       v, w, k_c, k_m: road parameters
        %       postm: length in meters
        %       start_time, end_time: start and end time of the simulation
        %       qin_meas: column vector, float, upstream measuremnt data
        %       qout_meas: column vector, float, downstream measuremnt data
        %       T_in: column vector, float, the upstream time grid, same
        %           length as qin_meas
        %       T_out: column vector, float, downstream time grid
        %       ini_segments: column vector, float, cumulative space grid
        %       rho_ini: column vector, float, initial density 
        %       e_max: the measurement error from the data
        function self = setIneqConstraints(...
                v, w, k_c, k_m, postm,...
                start_time, end_time,...
                qin_meas, qout_meas,...
                T_in, T_out,...
                ini_segments, rho_ini, e_max)            
            
            self.n_max_us = size(qin_meas,1);   % boundary conditions
            self.n_max_ds = size(qout_meas,1);
            
            self.b_max = size(rho_ini,1);   % initial conditions
            self.m_max = 0;                 % disable internal condition
            self.u_max = 0;                 % disable density condition
            
            self.v = v;
            self.w = w;
            self.k_c = k_c;
            self.k_m = k_m;
            
            self.start_time = start_time;
            self.end_time = end_time;
            
            self.xi = 0;  % relative postm for each link
            self.chi = postm;
            
            % if still using uniform discretization
            self.T_us = T_in;
            self.T_us_cum = [0; cumsum(T_in)];
            self.T_ds = T_out;
            self.T_ds_cum = [0; cumsum(T_out)];
            
            
            % set boundary condition
            self.qin_meas = qin_meas;
            self.qout_meas = qout_meas;
            
            % disable internal condition
            self.x_min = [];
            self.x_max = [];
            self.t_min = [];
            self.t_max = [];
            self.v_meas = [];
            self.r_meas = [];
            
            %disable density condition
            self.x_min_u = [];
            self.x_max_u = [];
            self.t_u = [];
            self.dens_meas = [];
            
            self.e_max = e_max;   
            
            % set initial condition
            if ~iscolumn(ini_segments)
                ini_segments = ini_segments';
            end
            self.ini_segments = ini_segments;
            
            if (length(ini_segments)~=length(rho_ini)+1 && ~isempty(rho_ini))
                error('Check ini_segments and rho_ini dimensions');
            end
            
            self.X = self.ini_segments(2:length(self.ini_segments)) - self.ini_segments(1:length(self.ini_segments)-1);
            
            if (isempty(rho_ini) || ~all(~isnan(rho_ini)))
                sprintf('Warning: Not defining a complete set of initial value conditions could significantly increase the computational time!\n');
            end
            
            if ~iscolumn(rho_ini)
                rho_ini = rho_ini';
            end
            self.rho_ini = rho_ini;
            
            if (~isempty(rho_ini))
                block_ini = double(rho_ini > k_c); 
                block_ini(isnan(rho_ini)) = NaN;
                self.indiced_rho_ini = self.groupSameElement(block_ini,inf);
            else
                self.indiced_rho_ini = [];
            end
            
            self.num_aux = 0;
            
            self.size_row = self.n_max_us + self.n_max_ds +...
                          self.b_max + 2*self.m_max + 2*self.u_max + 1; 
            %This row size if defined as a temporary one to compute the number of needed binary variables
            [self.nb_min, self.nb_max, self.nb_min_u, self.nb_max_u] = self.getBinaryvar;
            
            self.size_row = self.n_max_us + self.n_max_ds + self.b_max + 2*self.m_max + 2*self.u_max +...
                sum(self.nb_min) + sum(self.nb_max) + sum(self.nb_min_u) + sum(self.nb_max_u) + sum(self.num_aux) + 1;  
            self.list = [];
            
        end
        
        %==========================================================================
        %Method to substract to arrays with a "null" condition
        
        function [output] = substractArray(~,value1,value2)
            if ( (~isempty(value1))&& (~isempty(value2)) )
                output = value1-value2;
            else
                output = []; %return "null"
            end
        end
        
        %==========================================================================
        %Equations of function Mgamma, explicit solution of upstream condition
        
        function[array] = mgamma(self,n,t,x)
            
            array = zeros(1,self.size_row);  %Initialize the array
            
            if((self.T_us_cum(n) + (x - self.xi)/self.v <= t) &&...
                    (self.T_us_cum(n+1) + (x-self.xi)/self.v >= t ) && (n<= self.n_max_us))
                % in characteristic domain
                array(1,1:n-1) = self.T_us(1:n-1);
                array(1,n) = t - (x-self.xi) / self.v - self.T_us_cum(n);
                return
                
            elseif((self.T_us_cum(n+1) + (x-self.xi)/self.v < t) && (n <= self.n_max_us))
                % in Fan domain
                array(1,1:n) = self.T_us(1:n);
                array(1,self.size_row) = -self.k_c*self.v * (t - self.T_us_cum(n+1) - (x-self.xi)/self.v);
                return
                
            else
                array = [];     %return a "null"
                return
            end
            
        end
        
        %==========================================================================
        %Equations of function Mbeta, explicit solution of downstream condition
        
        function [array] = mbeta(self,n,t,x)
            array = zeros(1,self.size_row);
            
            if( (self.T_ds_cum(n) + (x - self.chi)/self.w <= t) &&...
                    (self.T_ds_cum(n+1) + (x-self.chi)/self.w >= t )&& (n <= self.n_max_ds))
                % in characteristic domain
                array(1,self.n_max_us + self.n_max_ds +1:...
                    self.n_max_us + self.n_max_ds + self.b_max) = -self.X; %initial number of vehicles
                array(1,self.n_max_us +1:self.n_max_us + n-1) = self.T_ds(1:n-1);
                array(1,self.n_max_us + n) = t - (x-self.chi)/self.w - self.T_ds_cum(n);
                array(1,self.size_row) = self.k_m*(x-self.chi);
                return
                
            elseif( (self.T_ds_cum(n+1) + (x-self.chi)/self.w < t) && (n <= self.n_max_ds))
                % in Fan domain
                array(1,self.n_max_us + self.n_max_ds +1:...
                    self.n_max_us + self.n_max_ds + self.b_max) = -self.X; %initial number of vehicles
                array(1,self.n_max_us +1:self.n_max_us + n) = self.T_ds(1:n);
                array(1,self.size_row) = -self.k_c*self.v * (t - self.T_ds_cum(n+1) - (x-self.chi)/(self.v));
                return
            else
                array = []; %return a "null"
                return
            end
        end
        
        %==========================================================================
        %Equations of function Mmu, explicit solutions of internal conditions
        
        function [array] = mmu(self, m, t, x)
            array = zeros(1,self.size_row);
            
            if( (x >= self.x_min(m) + self.v_meas(m)*(t-self.t_min(m))) &&...
                ( x >= self.x_max(m) + self.v*(t-self.t_max(m))) && (x <= self.x_min(m) + self.v*(t-self.t_min(m))) )
                array(1,self.n_max_us + self.n_max_ds + self.b_max + m) = 1;
                array(1,self.n_max_us + self.n_max_ds + self.b_max + self.m_max + m) =...
                    t - (x - self.x_min(m) - self.v_meas(m)*(t - self.t_min(m))) / (self.v - self.v_meas(m)) - self.t_min(m) ;
                return
            end
            
            if ( (x<= self.x_min(m) + self.v_meas(m)*(t-self.t_min(m))) &&...
                 ( x <= self.x_max(m) + self.w*(t-self.t_max(m))) && (x >= self.x_min(m) + self.w*(t-self.t_min(m))) )
                array(1,self.n_max_us + self.n_max_ds + self.b_max + m) = 1;
                array(1,self.n_max_us + self.n_max_ds + self.b_max + self.m_max + m) =...
                    t - (x - self.x_min(m) - self.v_meas(m)*(t-self.t_min(m)))/(self.w - self.v_meas(m)) - self.t_min(m);
                array(1,self.size_row) = -self.k_c*(self.v-self.w)*(x - self.x_min(m) - self.v_meas(m)*(t-self.t_min(m)))/(self.w - self.v_meas(m)) ;
                return
            end
            
            if( (x < self.x_max(m) + self.v*(t - self.t_max(m))) && ( x > self.x_max(m) + self.w*(t - self.t_max(m))) )
                array(1,self.n_max_us + self.n_max_ds + self.b_max + m) = 1;
                array(1, self.n_max_us + self.n_max_ds + self.b_max + self.m_max + m) = self.t_max(m) - self.t_min(m);
                array(1, self.size_row) = -(t-self.t_max(m))*self.k_c*(self.v - (x-self.x_max(m))/(t-self.t_max(m)+0.000001));
                return
            else
                array = []; %Return a "NULL"
                return
            end
        end
        
        %==========================================================================
        %Equations of function MTau1, explicit solutions of initial solutions applying when p<pc
        
        function [array] = mtau1 (self,b,t,x)
            array = zeros(1,self.size_row);
            
            if ( (self.xi + sum(self.X(1:b-1)) + t*self.v <= x) && (x <=self.xi+ sum(self.X(1:b)) + t*self.v) && ( b <= self.b_max))
                array(1,self.n_max_us + self.n_max_ds + 1 : self.n_max_us + self.n_max_ds + b-1) = -self.X(1:b-1);
                
                array(1,self.n_max_us + self.n_max_ds + b) = (t*self.v - x + sum(self.X(1:b-1)) + self.xi);
                return
            end
            
            if ( (self.xi+ sum(self.X(1:b-1)) + t*self.w <= x) && (x < self.xi+ sum(self.X(1:b-1)) + t*self.v) && (b <= self.b_max))
                array(1,self.n_max_us + self.n_max_ds + 1: self.n_max_us + self.n_max_ds + b-1) = -self.X(1:b-1);
                
                array(1,self.size_row) = -self.k_c*(t*self.v - x + sum(self.X(1:b-1)) + self.xi);
                return
            else
                array = []; %Return a "NULL"
                return
            end
        end
        
        %==========================================================================
        %Equations of function MTau2, explicit solutions of initial solutions applying when p>pc
        
        function [array] = mtau2 (self,b,t,x)
            array = zeros(1,self.size_row);
            if ( (self.xi + sum(self.X(1:b-1)) + t*self.w <= x) &&...
                    (x <= self.xi + sum(self.X(1:b)) + t*self.w) &&...
                    ( b <= self.b_max))
                
                array(1,self.n_max_us + self.n_max_ds + 1 : self.n_max_us + self.n_max_ds + b - 1 ) = -self.X(1:b-1);
                array(1,self.n_max_us + self.n_max_ds + b) = (t*self.w - x + sum(self.X(1:b-1)) + self.xi);
                array(1,self.size_row) = self.k_m*t*self.w;
                return
            end
            
            if ( (self.xi + sum(self.X(1:b)) + t*self.w < x) &&...
                    (x <= self.xi + sum(self.X(1:b)) + t*self.v) &&...
                    (b <= self.b_max))
                
                array(1,self.n_max_us + self.n_max_ds + 1:self.n_max_us + self.n_max_ds+b-1) = -self.X(1:b-1);
                array(1,self.n_max_us + self.n_max_ds + b) = -self.X(b);
                array(1,self.size_row) = -(self.k_c*(t*self.w - x + sum(self.X(1:b)) + self.xi) - self.k_m*t*self.w);
                return
            else
                array = []; %Return a "NULL" value
                return
            end
        end
        
        
        %==========================================================================        
        %Equations of function Mups1, explicit solutions of density constraints applying when p<pc        
        
        function [array] = mups1 (self,u,t,x)
            
            array = zeros(1,self.size_row);
            
            if ( (self.x_min_u(u) + (t-self.t_u(u))*self.v <= x) &&...
                 (x <=self.x_max_u(u) + (t-self.t_u(u))*self.v) &&...
                 (t>=self.t_u(u)) && ( u <= self.u_max))
                
                array(1,self.n_max_us + self.n_max_ds + self.b_max + 2*self.m_max + u ) = 1;
                array(1,self.n_max_us + self.n_max_ds + self.b_max + 2*self.m_max + self.u_max + u) =...
                    (t-self.t_u(u))*self.v - x + self.x_min_u(u);
                return
            end
            
            if ( (self.x_min_u(u) + (t-self.t_u(u))*self.w <= x) &&...
                 (x <= self.x_min_u(u) + (t-self.t_u(u))*self.v) &&...
                 (t>=self.t_u(u)) && (u <= self.u_max))
                array(1,self.n_max_us + self.n_max_ds + self.b_max + 2*self.m_max + u ) = 1;
                array(1,self.size_row) = -self.k_c*((t-self.t_u(u))*self.v - x + self.x_min_u(u));
                return
            else
                
                array = []; %Return a "NULL"
                return
                
            end
            
        end
      
        %==========================================================================        
        %Equations of function Mups2, explicit solutions of density constraints applying when p>pc        
        
        function [array] = mups2 (self,u,t,x)
            
            array = zeros(1,self.size_row);
            
            if ( (self.x_min_u(u) + (t-self.t_u(u))*self.w <= x) &&...
                 (x <=self.x_max_u(u) + (t-self.t_u(u))*self.w) &&...
                 (t>=self.t_u(u)) && ( u <= self.u_max))
                
                array(1,self.n_max_us + self.n_max_ds + self.b_max + 2*self.m_max + u ) = 1;
                array(1,self.n_max_us + self.n_max_ds + self.b_max + 2*self.m_max + self.u_max + u) =...
                        (t-self.t_u(u))*self.w - x + self.x_min_u(u);
                array(1,self.size_row) = self.k_m*(t-self.t_u(u))*self.w;
                return
            end
            
            if ( (self.x_max_u(u) + (t-self.t_u(u))*self.w <= x) &&...
                 (x <= self.x_max_u(u) + (t-self.t_u(u))*self.v) &&...
                 (t>=self.t_u(u)) && (u <= self.u_max))
                
                array(1,self.n_max_us + self.n_max_ds + self.b_max + 2*self.m_max + u ) = 1;
                array(1,self.size_row) = -(self.k_c*((t-self.t_u(u))*self.w - x +...
                                self.x_max_u(u)) - self.k_m*(t-self.t_u(u))*self.w);
                return
            else
                array = []; %Return a "NULL"
                return
            end
            
        end
        
        %==========================================================================
        %Equations that define the upstream condition
        
        function[array] = gamma(self,t,x)
            
            array = zeros(1,self.size_row);  %Initialize the array
            
            if t < self.start_time || t > self.end_time
                array = [];
                return
            end
            
            n = sum(self.T_us_cum <= t);  % in step n interval [ , )
            if t == self.end_time     % the final time point
                n = n-1;    % in the last step interval
            end
            
            
            if( (abs(x - self.xi)<0.00001) && (self.T_us_cum(n)<=t) &&...
                (t<=self.T_us_cum(n+1)) && (n <= self.n_max_us) )
                array(1,1:n-1) = self.T_us(1:n-1);
                array(1,n) = t-self.T_us_cum(n);
                return
            else
                array = []; %return a "NULL"
                return
            end
        end
        
        %==========================================================================
        %Equations that define the downstream condition
        
        function [array] = beta(self,t,x)
            array = zeros(1,self.size_row);
            
            if t < self.start_time || t > self.end_time
                array = [];
                return
            end
            
            n = sum(self.T_ds_cum <= t);  % in step n interval [ , )
            if t == self.end_time     % the final time point
                n = n-1;    % in the last step interval
            end
            
            if( (abs(x-self.chi)<0.01) && (self.T_ds_cum(n) <= t) &&...
                (t<= self.T_ds_cum(n+1)) && (n <= self.n_max_ds) )
                array(1,self.n_max_us + self.n_max_ds +1:...
                      self.n_max_us + self.n_max_ds + self.b_max) = -self.X;   %Initial number of vehicles
                array(1,self.n_max_us +1:self.n_max_us + n-1) = self.T_ds(1:n-1);
                array(1,self.n_max_us + n ) = t-self.T_ds_cum(n);
                return
            else
                array = [];
                return
            end
        end
        
        %==========================================================================
        %Equations that define the internal condition
        
        function [array] = mu(self, m ,t ,x)
            
            array = zeros(1,self.size_row);
            if ( (abs(x - (self.v_meas(m)*(t-self.t_min(m)) + self.x_min(m))) < 0.01) &&...
                    (t>= self.t_min(m)) && (t<= self.t_max(m)) )
                array(1, self.n_max_us + self.n_max_ds + self.b_max + m) = 1;
                array(1, self.n_max_us + self.n_max_ds + self.b_max + self.m_max + m ) = t - self.t_min(m);
                return
            else
                array = [];
                return
            end
        end
        
        %==========================================================================
        %Equations that define the density condition       

        %new variables x_min_u, x_max_u, t_u, u_max, m_max_vm0
      
        function [array] = ups(self, u ,t ,x)
            
            if ( (abs(t-self.t_u(u)) < 0.0001) && (x>= self.x_min_u(u)) && (x<= self.x_max_u(u)) )
                array = zeros(1,self.size_row);
                array(1, self.n_max_us + self.n_max_ds + self.b_max + 2*self.m_max + u) = 1;
                array(1, self.n_max_us + self.n_max_ds + self.b_max + 2*self.m_max + self.u_max + u ) = -(x - self.x_min_u(u+1));
                return
            else
                
                array = [];
                return
                
            end
            
        end
      
        %==========================================================================
        %Equations that define the initial condition
        
        function[array] = tau(self, b, t, x)
            if( (abs(t-self.t0)<0.01) &&  ( sum(self.X(1:b)) <= x) && ( x <=sum(self.X(1:b+1))) && (b <= self.b_max) )
                array = zeros(1,self.size_row);
                array(1, 2*(self.n_max+1)+1:2*(self.n_max+1)+b) = -self.X(1:b);
                array(1,2*(self.n_max+1) +  b +1) = -(x-sum(self.X(1:b)));
                return
            else
                array = []; %Return a "NULL"
                return
            end
        end
        
        %==========================================================================
        %Function to create the model constraints
        % output: the model constraints matrix
        
        function  [list] = setModelMatrix(self)
            
            % first allocate memory for list
            
            list = zeros(100000,self.size_row);
            % initialize rows counts
            rows = 0; 
            
            %==============================================
            % Solution associated with n-th upstream boundary conditions
            for n=1:self.n_max_us
                
                % at downstream points including the begining and ending
                for p=1:self.n_max_ds+1
                    
                    if self.T_us_cum(n) + (self.chi-self.xi)/self.v <= self.T_ds_cum(p) &&...
                       self.T_us_cum(n+1) + (self.chi-self.xi)/self.v >= self.T_ds_cum(p)
                        % point (self.T_ds_cum(p), self.chi)
                        array = self.mgamma(n,self.T_ds_cum(p),self.chi);
                        array2 = self.substractArray(array, self.beta(self.T_ds_cum(p), self.chi));
                        if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                            rows = rows+1;
                            list(rows,:) = array2;
                        end
                    end
                    
                end
                
                % freeflow speed intersection at downstream point
                % point (self.T_us_cum(n) + (self.chi-self.xi)/self.v, self.chi) 
                % where solution >= value condition 
                array = self.mgamma(n, self.T_us_cum(n) + (self.chi-self.xi)/self.v, self.chi);
                array2 = self.substractArray(array, self.beta(self.T_us_cum(n) + (self.chi-self.xi)/self.v, self.chi));
                if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                    rows = rows+1;
                    list(rows,:) = array2;
                end
                
                % There is no need to pose the inequality for upstream
                % points or initial time points
                
                % internal condition points
                for m=1:self.m_max
                    
                    array = self.mgamma(n,self.t_min(m), self.x_min(m));
                    array2 = self.substractArray(array, self.mu(m, self.t_min(m), self.x_min(m)));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        
                        list(rows,:) = array2;
                    end
                    
                    array = self.mgamma(n,self.t_max(m), self.x_max(m));
                    array2 = self.substractArray(array, self.mu(m, self.t_max(m), self.x_max(m)));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.T_us_cum(n) * self.v - self.v_meas(m) * self.t_min(m) +...
                              self.x_min(m) - self.xi) / (self.v - self.v_meas(m));
                    x_temp = self.x_min(m) + self.v_meas(m)*(t_temp - self.t_min(m));
                    array = self.mgamma(n,t_temp, x_temp);
                    array2 = self.substractArray(array, self.mu(m, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
          
                end
                
                % density condition points
                for u=1:self.u_max
                    
                    array = self.mgamma(n,self.t_u(u), self.x_min_u(u));
                    array2 = self.substractArray(array, self.ups(u, self.t_u(u), self.x_min_u(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    
                    array = self.mgamma(n,self.t_u(u), self.x_max_u(u));
                    array2 = self.substractArray(array, self.ups(u, self.t_u(u), self.x_max_u(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.mgamma(n,self.t_u(u), self.xi + self.v*(self.t_u(u)-self.T_us_cum(n)));
                    array2 = self.substractArray(array, self.ups(u, self.t_u(u), self.xi +...
                        self.v*(self.t_u(u)-self.T_us_cum(n))));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
     
                end
                
            end
            
            
            %==============================================
            % Solution associated with n-th downstram boundary conditions
            for n=1:self.n_max_ds
                
                % at upstream points
                for p=1:self.n_max_us+1
                    
                    if self.T_ds_cum(n) + (self.xi-self.chi)/self.w <= self.T_us_cum(p) &&...
                            self.T_ds_cum(n+1) + (self.xi-self.chi)/self.w >= self.T_us_cum(p)
                        % point (self.T_us_cum(p), self.xi)
                        array = self.mbeta(n,self.T_us_cum(p), self.xi);
                        array2 = self.substractArray(array, self.gamma(self.T_us_cum(p), self.xi));
                        if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                            rows = rows+1;
                            
                            list(rows,:) = array2;
                        end
                    end
                    
                end
                
                % intersection point at upstream (self.T_ds_cum(n) + (self.xi-self.chi)/self.w,self.xi)
                array = self.mbeta(n,self.T_ds_cum(n) + (self.xi-self.chi)/self.w, self.xi);
                array2 = self.substractArray(array, self.gamma(self.T_ds_cum(n) + (self.xi-self.chi)/self.w, self.xi));
                if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                    rows = rows+1;
                    
                    list(rows,:) = array2;
                end
                
                % There is no need to pose the inequality for upstream
                % points or initial time points
                
                % internal condition points
                for m=1:self.m_max
                    array = self.mbeta(n,self.t_min(m), self.x_min(m));
                    array2 = self.substractArray(array, self.mu(m, self.t_min(m), self.x_min(m)));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.mbeta(n,self.t_max(m), self.x_max(m));
                    array2 = self.substractArray(array, self.mu(m, self.t_max(m), self.x_max(m)));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.T_ds_cum(n) * self.w - self.v_meas(m) * self.t_min(m) +...
                        self.x_min(m) - self.chi) / (self.w - self.v_meas(m));
                    x_temp = self.x_min(m) + self.v_meas(m)*(t_temp - self.t_min(m));
                    array = self.mbeta(n,t_temp, x_temp);
                    array2 = self.substractArray(array, self.mu(m, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                end
                
                % density condition points
                for u=1:self.u_max
                    
                    array = self.mbeta(n,self.t_u(u), self.x_min_u(u));
                    array2 = self.substractArray(array, self.ups(u, self.t_u(u), self.x_min_u(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.mbeta(n,self.t_u(u), self.x_max_u(u));
                    array2 = self.substractArray(array, self.ups(u, self.t_u(u), self.x_max_u(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.mbeta(n,self.t_u(u), self.chi + self.w*(self.t_u(u)-self.T_ds_cum(n)));
                    array2 = self.substractArray(array, self.ups(u, self.t_u(u), self.chi +...
                        self.w*(self.t_u(u)-self.T_ds_cum(n))));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                end
                
            end
            
            
            %==============================================
            % Solution associated with initial conditions
            for k=1:self.b_max
                
                % at upstream points
                for p=1:self.n_max_us+1
                    
                    %TAU1 rho < rho_c
                    
                    % upstream time grid points
                    array = self.mtau1(k,self.T_us_cum(p),self.xi);
                    array2 = self.substractArray(array,self.gamma(self.T_us_cum(p),self.xi));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    
                    %TAU 2
                    array = self.mtau2(k,self.T_us_cum(p),self.xi);
                    array2 = self.substractArray(array,self.gamma(self.T_us_cum(p),self.xi));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                end
                
                % at upstream intersection points
                % tau 1
                array = self.mtau1(k,self.start_time + (-sum(self.X(1:k)))/self.w,self.xi);
                array2 = self.substractArray(array,self.gamma(self.start_time +...
                    (-sum(self.X(1:k)))/self.w,self.xi));
                if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                    rows = rows+1;
                    list(rows,:) = array2;
                end
                
                % tau 2
                array = self.mtau2(k,self.start_time + (-sum(self.X(1:k)))/self.w,self.xi);
                array2 = self.substractArray(array,self.gamma(self.start_time +...
                    (-sum(self.X(1:k)))/self.w,self.xi));
                if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                    rows = rows+1;
                    list(rows,:) = array2;
                end
                
                
                % at downstream points
                for p=1:self.n_max_ds+1
                   
                    % Tau1
                    array = self.mtau1(k,self.T_ds_cum(p),self.chi);
                    array2 = self.substractArray(array,self.beta(self.T_ds_cum(p),self.chi));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                   
                    %tau2
                    array = self.mtau2(k,self.T_ds_cum(p),self.chi);
                    array2 = self.substractArray(array,self.beta(self.T_ds_cum(p),self.chi));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                end
                
                % intersection point at downstream
                array = self.mtau1(k,self.start_time + (self.chi-(sum(self.X(1:k-1))+self.xi))/self.v,self.chi);
                array2 = self.substractArray(array,self.beta(self.start_time +...
                    (self.chi-(sum(self.X(1:k-1))+self.xi))/self.v,self.chi));
                if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                    rows = rows+1;
                    list(rows,:) = array2;
                end
                
                array = self.mtau2(k,self.start_time + (self.chi-(sum(self.X(1:k-1))+self.xi))/self.v,self.chi);
                array2 = self.substractArray(array,self.beta(self.start_time +...
                    (self.chi-(sum(self.X(1:k-1))+self.xi))/self.v,self.chi));
                if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                    rows = rows+1;
                    list(rows,:) = array2;
                end
                
                % for internal points
                for p=1:self.m_max
                    
                    %TAU 1
                    array = self.mtau1(k,self.t_min(p), self.x_min(p));
                    array2 = self.substractArray(array, self.mu(p, self.t_min(p), self.x_min(p)));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.mtau1(k,self.t_max(p), self.x_max(p));
                    array2 = self.substractArray(array, self.mu(p, self.t_max(p), self.x_max(p)));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.x_min(p) - (self.xi + sum(self.X(1:k))) - self.v_meas(p) * self.t_min(p)) ...
                        / (self.v - self.v_meas(p));
                    x_temp = self.v_meas(p) * (t_temp - self.t_min(p)) + self.x_min(p);
                    array = self.mtau1(k,t_temp, x_temp);
                    array2 = self.substractArray(array, self.mu(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.x_min(p) - (self.xi + sum(self.X(1:k-1))) - self.v_meas(p) * self.t_min(p))...
                        / (self.w - self.v_meas(p));
                    x_temp = self.v_meas(p) * (t_temp - self.t_min(p)) + self.x_min(p);
                    array = self.mtau1(k,t_temp, x_temp);
                    array2 = self.substractArray(array, self.mu(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.x_min(p) - (self.xi + sum(self.X(1:k-1))) - self.v_meas(p) * self.t_min(p))...
                        / (self.v - self.v_meas(p));
                    x_temp = self.v_meas(p) * (t_temp - self.t_min(p)) + self.x_min(p);
                    array = self.mtau1(k,t_temp, x_temp);
                    array2 = self.substractArray(array, self.mu(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.x_min(p) - (self.xi + sum(self.X(1:k))) - self.v_meas(p) * self.t_min(p))...
                        / (self.w - self.v_meas(p));
                    x_temp = self.v_meas(p) * (t_temp - self.t_min(p)) + self.x_min(p);
                    array = self.mtau1(k,t_temp, x_temp);
                    array2 = self.substractArray(array, self.mu(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    
                    %TAU 2
                    array = self.mtau2(k,self.t_min(p), self.x_min(p));
                    array2 = self.substractArray(array, self.mu(p, self.t_min(p), self.x_min(p)));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.mtau2(k,self.t_max(p), self.x_max(p));
                    array2 = self.substractArray(array, self.mu(p, self.t_max(p), self.x_max(p)));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.x_min(p) - (self.xi + sum(self.X(1:k))) - self.v_meas(p) * self.t_min(p))...
                        / (self.v - self.v_meas(p));
                    x_temp = self.v_meas(p) * (t_temp - self.t_min(p)) + self.x_min(p);
                    array = self.mtau2(k,t_temp, x_temp);
                    array2 = self.substractArray(array, self.mu(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.x_min(p) - (self.xi + sum(self.X(1:k-1))) - self.v_meas(p) * self.t_min(p))...
                        / (self.w - self.v_meas(p));
                    x_temp = self.v_meas(p) * (t_temp - self.t_min(p)) + self.x_min(p);
                    array = self.mtau2(k,t_temp, x_temp);
                    array2 = self.substractArray(array, self.mu(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.x_min(p) - (self.xi + sum(self.X(1:k-1))) - self.v_meas(p) * self.t_min(p))...
                        / (self.v - self.v_meas(p));
                    x_temp = self.v_meas(p) * (t_temp - self.t_min(p)) + self.x_min(p);
                    array = self.mtau2(k,t_temp, x_temp);
                    array2 = self.substractArray(array, self.mu(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.x_min(p) - (self.xi + sum(self.X(1:k))) - self.v_meas(p) * self.t_min(p))...
                        / (self.w - self.v_meas(p));
                    x_temp = self.v_meas(p) * (t_temp - self.t_min(p)) + self.x_min(p);
                    array = self.mtau2(k,t_temp, x_temp);
                    array2 = self.substractArray(array, self.mu(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                end
                
                % for density points
                for u=1:self.u_max
                    
                    %TAU 1
                    array = self.mtau1(k,self.t_u(u), self.x_min_u(u));
                    array2 = self.substractArray(array, self.ups(u, self.t_u(u), self.x_min_u(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.mtau1(k,self.t_u(u), self.x_max_u(u));
                    array2 = self.substractArray(array, self.ups(u, self.t_u(u), self.x_max_u(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.mtau1(k,self.t_u(u), self.xi + sum(self.X(1:k-1)) +...
                        (self.t_u(u)-self.start_time)*self.v);
                    array2 = self.substractArray(array, self.ups(u, self.t_u(u), self.xi + sum(self.X(1:k-1)) +...
                        (self.t_u(u)-self.start_time)*self.v));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.mtau1(k,self.t_u(u), self.xi + sum(self.X(1:k)) +...
                        (self.t_u(u)-self.start_time)*self.v);
                    array2 = self.substractArray(array, self.ups(u, self.t_u(u), self.xi + sum(self.X(1:k)) +...
                        (self.t_u(u)-self.start_time)*self.v));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.mtau1(k,self.t_u(u), self.xi + sum(self.X(1:k-1)) +...
                        (self.t_u(u)-self.start_time)*self.w);
                    array2 = self.substractArray(array, self.ups(u, self.t_u(u), self.xi + sum(self.X(1:k-1)) +...
                        (self.t_u(u)-self.start_time)*self.w));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.mtau1(k,self.t_u(u), self.xi + sum(self.X(1:k)) +...
                        (self.t_u(u)-self.start_time)*self.w);
                    array2 = self.substractArray(array, self.ups(u, self.t_u(u), self.xi + sum(self.X(1:k)) +...
                        (self.t_u(u)-self.start_time)*self.w));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    
                    %TAU 2
                    array = self.mtau2(k,self.t_u(u), self.x_min_u(u));
                    array2 = self.substractArray(array, self.ups(u, self.t_u(u), self.x_min_u(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.mtau2(k,self.t_u(u), self.x_max_u(u));
                    array2 = self.substractArray(array, self.ups(u, self.t_u(u), self.x_max_u(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.mtau2(k,self.t_u(u), self.xi + sum(self.X(1:k-1)) +...
                        (self.t_u(u)-self.start_time)*self.v);
                    array2 = self.substractArray(array, self.ups(u, self.t_u(u), self.xi + sum(self.X(1:k-1)) +...
                        (self.t_u(u)-self.start_time)*self.v));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.mtau2(k,self.t_u(u), self.xi + sum(self.X(1:k)) +...
                        (self.t_u(u)-self.start_time)*self.v);
                    array2 = self.substractArray(array, self.ups(u, self.t_u(u), self.xi + sum(self.X(1:k)) +...
                        (self.t_u(u)-self.start_time)*self.v));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.mtau2(k,self.t_u(u), self.xi + sum(self.X(1:k-1)) +...
                        (self.t_u(u)-self.start_time)*self.w);
                    array2 = self.substractArray(array, self.ups(u, self.t_u(u), self.xi + sum(self.X(1:k-1)) +...
                        (self.t_u(u)-self.start_time)*self.w));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.mtau2(k,self.t_u(u), self.xi + sum(self.X(1:k)) +...
                        (self.t_u(u)-self.start_time)*self.w);
                    array2 = self.substractArray(array, self.ups(u, self.t_u(u), self.xi + sum(self.X(1:k)) +...
                        (self.t_u(u)-self.start_time)*self.w));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                end
                
            end
            
            
            %==============================================
            % Solution associated with internal conditions
            for m=1:self.m_max
                
                % upstream points
                for p=1:self.n_max_us+1
                    
                    array = self.mmu(m,self.T_us_cum(p), self.xi);
                    array2 = self.substractArray(array, self.gamma(self.T_us_cum(p), self.xi));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.xi - self.x_min(m) + self.w*self.t_min(m)) / (self.w);
                    x_temp = self.xi;
                    array = self.mmu(m,t_temp, x_temp);
                    array2 = self.substractArray(array, self.gamma(t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.xi - self.x_max(m) + self.w*self.t_max(m)) / (self.w);
                    x_temp = self.xi;
                    array = self.mmu(m,t_temp, x_temp);
                    array2 = self.substractArray(array, self.gamma(t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                end
                
                
                % downstream points
                for p=1:self.n_max_ds+1
                    
                    array = self.mmu(m,self.T_ds_cum(p), self.chi);
                    array2 = self.substractArray(array, self.beta(self.T_ds_cum(p), self.chi));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.chi - self.x_min(m) + self.v*self.t_min(m)) / (self.v);
                    x_temp = self.chi;
                    array = self.mmu(m,t_temp, x_temp);
                    array2 = self.substractArray(array, self.beta(t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.chi - self.x_max(m) + self.v*self.t_max(m)) / (self.v);
                    x_temp = self.chi;
                    array = self.mmu(m,t_temp, x_temp);
                    array2 = self.substractArray(array, self.beta(t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                end
                
                
                % other internal points
                for p=1:self.m_max
                    
                    array = self.mmu(m,self.t_min(p), self.x_min(p));
                    array2 = self.substractArray(array, self.mu(p,self.t_min(p), self.x_min(p)));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.mmu(m,self.t_max(p), self.x_max(p));
                    array2 = self.substractArray(array, self.mu(p,self.t_max(p), self.x_max(p)));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.x_min(m) - self.x_min(p) + self.v_meas(p) * self.t_min(p) -...
                        self.v_meas(m) * self.t_min(m)) / (self.v_meas(p) - self.v_meas(m));
                    x_temp = self.v_meas(p) * (t_temp - self.t_min(p)) + self.x_min(p);
                    array = self.mmu(m,t_temp, x_temp);
                    array2 = self.substractArray(array, self.mu(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    
                    t_temp = (self.x_max(m) - self.x_min(p) + self.v_meas(p) * self.t_min(p) -...
                        self.v * self.t_max(m)) / (self.v_meas(p) - self.v);
                    x_temp = self.v_meas(p) * (t_temp - self.t_min(p)) + self.x_min(p);
                    array = self.mmu(m,t_temp, x_temp);
                    array2 = self.substractArray(array, self.mu(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    
                    t_temp = (self.x_min(m) - self.x_min(p) + self.v_meas(p) * self.t_min(p) -...
                        self.v * self.t_min(m)) / (self.v_meas(p) - self.v);
                    x_temp = self.v_meas(p) * (t_temp - self.t_min(p)) + self.x_min(p);
                    array = self.mmu(m,t_temp, x_temp);
                    array2 = self.substractArray(array, self.mu(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    
                    t_temp = (self.x_max(m) - self.x_min(p) + self.v_meas(p) * self.t_min(p) -...
                        self.v * self.t_max(m)) / (self.v_meas(p) - self.w);
                    x_temp = self.v_meas(p) * (t_temp - self.t_min(p)) + self.x_min(p);
                    array = self.mmu(m,t_temp, x_temp);
                    array2 = self.substractArray(array, self.mu(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    
                    t_temp = (self.x_min(m) - self.x_min(p) + self.v_meas(p) * self.t_min(p) -...
                        self.v * self.t_min(m)) / (self.v_meas(p) - self.w);
                    x_temp = self.v_meas(p) * (t_temp - self.t_min(p)) + self.x_min(p);
                    array = self.mmu(m,t_temp, x_temp);
                    array2 = self.substractArray(array, self.mu(p, t_temp, x_temp));
                    if(~isempty(array2) && ~all(array2<=10e-11 & array2 >=-10e-11))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
      
                    
                end
                
                
                % density condition points
                for u=1:self.u_max
                    
                    array = self.mmu(m,self.t_u(u), self.x_min_u(u));
                    array2 = self.substractArray(array, self.ups(u,self.t_u(u), self.x_min_u(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.mmu(m,self.t_u(u), self.x_max_u(u));
                    array2 = self.substractArray(array, self.ups(u,self.t_u(u), self.x_max_u(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.mmu(m,self.t_u(u), self.x_min(m)+(self.t_u(u)-self.t_min(m))*self.w);
                    array2 = self.substractArray(array, self.ups(u,self.t_u(u),...
                                               self.x_min(m)+(self.t_u(u)-self.t_min(m))*self.w));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.mmu(m,self.t_u(u), self.x_max(m)+(self.t_u(u)-self.t_max(m))*self.w);
                    array2 = self.substractArray(array, self.ups(u,self.t_u(u),...
                        self.x_max(m)+(self.t_u(u)-self.t_max(m))*self.w));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.mmu(m,self.t_u(u), self.x_min(m)+(self.t_u(u)-self.t_min(m))*self.v);
                    array2 = self.substractArray(array, self.ups(u,self.t_u(u),...
                        self.x_min(m)+(self.t_u(u)-self.t_min(m))*self.v));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.mmu(m,self.t_u(u), self.x_max(m)+(self.t_u(u)-self.t_max(m))*self.v);
                    array2 = self.substractArray(array, self.ups(u,self.t_u(u),...
                        self.x_max(m)+(self.t_u(u)-self.t_max(m))*self.v));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                end
                
            end
            
            
            % Solution associated with density conditions
            for k=1:self.u_max

                % upstream points
                for p=1:self.n_max_us+1
                    
                    %UPS1
                    array = self.mups1(k,self.T_us_cum(p),self.xi);
                    array2 = self.substractArray(array,self.gamma(self.T_us_cum(p),self.xi));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    %UPS2
                    array = self.mups2(k,self.T_us_cum(p),self.xi);
                    array2 = self.substractArray(array,self.gamma(self.T_us_cum(p),self.xi));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                end
                
                % upstream intersection points
                array = self.mups1(k,self.t_u(k) + (self.xi-self.x_max_u(k))/self.w,self.xi);
                array2 = self.substractArray(array,self.gamma(self.t_u(k) +...
                        (self.xi-self.x_max_u(k))/self.w,self.xi));
                if(~isempty(array2))
                    rows = rows+1;
                    list(rows,:) = array2;
                end
                
                array = self.mups2(k,self.t_u(k) + (self.xi-self.x_max_u(k))/self.w,self.xi);
                array2 = self.substractArray(array,self.gamma(self.t_u(k) +...
                    (self.xi-self.x_max_u(k))/self.w,self.xi));
                if(~isempty(array2))
                    rows = rows+1;
                    list(rows,:) = array2;
                end
                
                
                % downstream points
                for p=1:self.n_max_ds+1
                    
                    %UPS1
                    array = self.mups1(k,self.T_ds_cum(p),self.chi);
                    array2 = self.substractArray(array,self.beta(self.T_ds_cum(p),self.chi));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    %USP2
                    array = self.mups2(k,self.T_ds_cum(p),self.chi);
                    array2 = self.substractArray(array,self.beta(self.T_ds_cum(p),self.chi));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                end
                
                % downstream intersection points
                array = self.mups1(k,self.t_u(k) + (self.chi-self.x_min_u(k))/self.v,self.chi);
                array2 = self.substractArray(array,self.beta(self.t_u(k) +...
                    (self.chi-self.x_min_u(k))/self.v,self.chi));
                if(~isempty(array2))
                    rows = rows+1;
                    list(rows,:) = array2;
                end
                
                array = self.mups2(k,self.t_u(k) + (self.chi-self.x_min_u(k))/self.v,self.chi);
                array2 = self.substractArray(array,self.beta(self.t_u(k) +...
                    (self.chi-self.x_min_u(k))/self.v,self.chi));
                if(~isempty(array2))
                    rows = rows+1;
                    list(rows,:) = array2;
                end
                
                
                % internal condition points
                for p=1:self.m_max
                    
                    %UPS 1
                    array = self.mups1(k,self.t_min(p), self.x_min(p));
                    array2 = self.substractArray(array, self.mu(p, self.t_min(p), self.x_min(p)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.mups1(k,self.t_max(p), self.x_max(p));
                    array2 = self.substractArray(array, self.mu(p, self.t_max(p), self.x_max(p)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.x_min(p) - self.x_max_u(k) - self.v_meas(p) * self.t_min(p) +...
                        self.v*self.t_u(k)) / (self.v - self.v_meas(p));
                    x_temp = self.v_meas(p) * (t_temp - self.t_min(p)) + self.x_min(p);
                    array = self.mups1(k,t_temp, x_temp);
                    array2 = self.substractArray(array, self.mu(p, t_temp, x_temp));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.x_min(p) - self.x_min_u(k) - self.v_meas(p) * self.t_min(p) +...
                        self.w*self.t_u(k)) / (self.w - self.v_meas(p));
                    x_temp = self.v_meas(p) * (t_temp - self.t_min(p)) + self.x_min(p);
                    array = self.mups1(k,t_temp, x_temp);
                    array2 = self.substractArray(array, self.mu(p, t_temp, x_temp));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.x_min(p) - self.x_min_u(k) - self.v_meas(p) * self.t_min(p) +...
                        self.v*self.t_u(k)) / (self.v - self.v_meas(p));
                    x_temp = self.v_meas(p) * (t_temp - self.t_min(p)) + self.x_min(p);
                    array = self.mups1(k,t_temp, x_temp);
                    array2 = self.substractArray(array, self.mu(p, t_temp, x_temp));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.x_min(p) - self.x_max_u(k) - self.v_meas(p) * self.t_min(p) +...
                        self.w*self.t_u(k)) / (self.w - self.v_meas(p));
                    x_temp = self.v_meas(p) * (t_temp - self.t_min(p)) + self.x_min(p);
                    array = self.mups1(k,t_temp, x_temp);
                    array2 = self.substractArray(array, self.mu(p, t_temp, x_temp));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    
                    %UPS 2
                    array = self.mups2(k,self.t_min(p), self.x_min(p));
                    array2 = self.substractArray(array, self.mu(p, self.t_min(p), self.x_min(p)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.mups2(k,self.t_max(p), self.x_max(p));
                    array2 = self.substractArray(array, self.mu(p, self.t_max(p), self.x_max(p)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.x_min(p) - self.x_max_u(k) - self.v_meas(p) * self.t_min(p) +...
                        self.v*self.t_u(k)) / (self.v - self.v_meas(p));
                    x_temp = self.v_meas(p) * (t_temp - self.t_min(p)) + self.x_min(p);
                    array = self.mups2(k,t_temp, x_temp);
                    array2 = self.substractArray(array, self.mu(p, t_temp, x_temp));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.x_min(p) - self.x_min_u(k) - self.v_meas(p) * self.t_min(p) +...
                        self.w*self.t_u(k)) / (self.w - self.v_meas(p));
                    x_temp = self.v_meas(p) * (t_temp - self.t_min(p)) + self.x_min(p);
                    array = self.mups2(k,t_temp, x_temp);
                    array2 = self.substractArray(array, self.mu(p, t_temp, x_temp));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.x_min(p) - self.x_min_u(k) - self.v_meas(p) * self.t_min(p) +...
                        self.v*self.t_u(k)) / (self.v - self.v_meas(p));
                    x_temp = self.v_meas(p) * (t_temp - self.t_min(p)) + self.x_min(p);
                    array = self.mups2(k,t_temp, x_temp);
                    array2 = self.substractArray(array, self.mu(p, t_temp, x_temp));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    t_temp = (self.x_min(p) - self.x_max_u(k) - self.v_meas(p) * self.t_min(p) +...
                        self.w*self.t_u(k)) / (self.w - self.v_meas(p));
                    x_temp = self.v_meas(p) * (t_temp - self.t_min(p)) + self.x_min(p);
                    array = self.mups2(k,t_temp, x_temp);
                    array2 = self.substractArray(array, self.mu(p, t_temp, x_temp));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                end
                
                
                % density condition points
                for u=1:self.u_max
                    
                    %UPS 1
                    array = self.mups1(k,self.t_u(u), self.x_min_u(u));
                    array2 = self.substractArray(array, self.ups(u, self.t_u(u), self.x_min_u(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.mups1(k,self.t_u(u), self.x_max_u(u));
                    array2 = self.substractArray(array, self.ups(u, self.t_u(u), self.x_max_u(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.mups1(k,self.t_u(u), self.x_min_u(k) + (self.t_u(u)-self.t_u(k))*self.v);
                    array2 = self.substractArray(array, self.ups(u, self.t_u(u), self.x_min_u(k) +...
                        (self.t_u(u)-self.t_u(k))*self.v));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.mups1(k,self.t_u(u), self.x_max_u(k) + (self.t_u(u)-self.t_u(k))*self.v);
                    array2 = self.substractArray(array, self.ups(u, self.t_u(u), self.x_max_u(k) +...
                        (self.t_u(u)-self.t_u(k))*self.v));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.mups1(k,self.t_u(u), self.x_min_u(k) + (self.t_u(u)-self.t_u(k))*self.w);
                    array2 = self.substractArray(array, self.ups(u, self.t_u(u), self.x_min_u(k) +...
                        (self.t_u(u)-self.t_u(k))*self.w));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.mups1(k,self.t_u(u), self.x_max_u(k) + (self.t_u(u)-self.t_u(k))*self.w);
                    array2 = self.substractArray(array, self.ups(u, self.t_u(u), self.x_max_u(k) +...
                        (self.t_u(u)-self.t_u(k))*self.w));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    % %UPS 2
                    array = self.mups2(k,self.t_u(u), self.x_min_u(u));
                    array2 = self.substractArray(array, self.ups(u, self.t_u(u), self.x_min_u(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.mups2(k,self.t_u(u), self.x_max_u(u));
                    array2 = self.substractArray(array, self.ups(u, self.t_u(u), self.x_max_u(u)));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.mups2(k,self.t_u(u), self.x_min_u(k) + (self.t_u(u)-self.t_u(k))*self.v);
                    array2 = self.substractArray(array, self.ups(u, self.t_u(u), self.x_min_u(k) +...
                        (self.t_u(u)-self.t_u(k))*self.v));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.mups2(k,self.t_u(u), self.x_max_u(k) + (self.t_u(u)-self.t_u(k))*self.v);
                    array2 = self.substractArray(array, self.ups(u, self.t_u(u), self.x_max_u(k) +...
                        (self.t_u(u)-self.t_u(k))*self.v));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.mups2(k,self.t_u(u), self.x_min_u(k) + (self.t_u(u)-self.t_u(k))*self.w);
                    array2 = self.substractArray(array, self.ups(u, self.t_u(u), self.x_min_u(k) +...
                        (self.t_u(u)-self.t_u(k))*self.w));
                    if(~isempty(array2))
                        rows = rows+1;
                        list(rows,:) = array2;
                    end
                    
                    array = self.mups2(k,self.t_u(u), self.x_max_u(k) + (self.t_u(u)-self.t_u(k))*self.w);
                    array2 = self.substractArray(array, self.ups(u, self.t_u(u), self.x_max_u(k) +...
                        (self.t_u(u)-self.t_u(k))*self.w));
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
        %Function to create the data constraints
        % output: the data constraints matrix
        
        function [list2] = setDataMatrix(self)
                        
            %Use sparse matrix to speed up computation
            list2 = zeros(10000,self.size_row);
            rows = 0;   %initialize row counts
            
            % upstream data constraints
            if (~isempty(self.qin_meas))
                
                for n=1:self.n_max_us
                    if ~isnan(self.qin_meas(n))
                        array = zeros(1,self.size_row);
                        array(1,n) = 1;
                        array(1,self.size_row) = self.qin_meas(n)*(1-self.e_max);
                        
                        rows = rows+1;
                        list2(rows,:)=array;
                    end
                end
                
                for n=1:self.n_max_us
                    if ~isnan(self.qin_meas(n))
                    array = zeros(1,self.size_row);
                    array(1,n) = -1;
                    array(1,self.size_row) = -self.qin_meas(n)*(1+self.e_max);
                    rows = rows+1;
                    list2(rows,:)=array;
                    end
                end
                
            end
            
            % downstream data constraints
            if (~isempty(self.qout_meas) )
                
                for n=1:self.n_max_ds
                    if ~isnan(self.qout_meas(n))
                    array = zeros(1,self.size_row);
                    array(1,self.n_max_us + n) = 1;
                    array(1,self.size_row) = self.qout_meas(n)*(1-self.e_max);
                    rows = rows+1;
                    list2(rows,:)=array;
                    end
                end
                
                for n=1:self.n_max_ds
                    if ~isnan(self.qout_meas(n))
                    array = zeros(1,self.size_row);
                    array(1,self.n_max_ds + n) = -1;
                    array(1,self.size_row) = -self.qout_meas(n)*(1+self.e_max);
                    rows = rows+1;
                    list2(rows,:)=array;
                    end
                end
                
            end
            
            % Initial data constraints
            if(~isempty(self.rho_ini) )
                
                for n=1:self.b_max
                    if ~isnan(self.rho_ini(n))
                        array = zeros(1,self.size_row);
                        array(1,self.n_max_us + self.n_max_ds + n) = 1;
                        array(1,self.size_row) = self.rho_ini(n)*(1-self.e_max);
                        rows = rows+1;
                        list2(rows,:)=array;
                    end
                end
                
                for n=1:self.b_max
                    if ~isnan(self.rho_ini(n))
                        array = zeros(1,self.size_row);
                        array(1,self.n_max_us + self.n_max_ds + n) = -1;
                        array(1,self.size_row) = -self.rho_ini(n)*(1+self.e_max);
                        rows = rows+1;
                        list2(rows,:)=array;
                    end
                end
            end
            
            % internal data constraints for rate of change r
            if(~isempty(self.r_meas))
                
                for m = 1:self.m_max
                    if ~isnan(self.r_meas(m))
                        array = zeros(1,self.size_row);
                        array(1,self.n_max_us + self.n_max_ds + self.b_max + self.m_max + m) = 1;
                        array(1,self.size_row) = self.r_meas(m);
                        rows = rows+1;
                        list2(rows,:)=array;
                    end
                end
                
                for m = 1:self.m_max
                    if ~isnan(self.r_meas(m))
                        array = zeros(1,self.size_row);
                        array(1,self.n_max_us + self.n_max_ds + self.b_max + self.m_max + m) = -1;
                        array(1,self.size_row) = - self.r_meas(m);
                        rows = rows+1;
                        list2(rows,:)=array;
                    end
                end
                
            end
            
            % density data constraints for density rho_u
            if(~isempty(self.dens_meas))
                
                for u = 1:self.u_max
                    if ~isnan(self.dens_meas(u))
                        array = zeros(1,self.size_row);
                        array(1,self.n_max_us + self.n_max_ds + self.b_max + 2*self.m_max + ...
                                self.u_max + u) = 1;
                        array(1,self.size_row) = self.dens_meas(u)*(1-self.e_max);
                        rows = rows+1;
                        list2(rows,:)=array;
                    end
                end
                
                for u = 1:self.u_max
                    if ~isnan(self.dens_meas(u))
                        array = zeros(1,self.size_row);
                        array(1,self.n_max_us + self.n_max_ds + self.b_max + 2*self.m_max + ...
                                self.u_max + u) = -1;
                        array(1,self.size_row) = -self.dens_meas(u)*(1+self.e_max);
                        rows = rows+1;
                        list2(rows,:)=array;
                    end
                end
                
            end
            
            
            % truncate matrix
            list2 = list2(1:rows,:);
            
            
        end
        
        %==========================================================================
        %Function to create the MILP constraints
        % The following functions need to be modified to enable different
        % discretization in up/downstream boundaries.
        
%         function [list3] = setMatrix3(si)
%             
%             list3 = zeros(0,0);
%             
%             %Binary variables useful counter as a reference
%             countB=0;
%             Cmax = 500000; 
%             
%             % nb: number of binary variables per point
%             nb = ceil(log2(2*(self.b_max+1)+2*(self.n_max+1)+(self.m_max+1)+(self.u_max+1)));
%             
%             % comb_mat: combination matrix for possible binary combinations
%             comb_mat = zeros(2^nb,nb);
%             
%             for i = 1:2^nb
%                 comb_mat(i,:) = double(dec2bin(i-1,nb)=='1');
%             end
%             
%             for k = 0:self.m_max
%                 
%                 %x_min and t_min point------------------
%                 if (k==0 || (self.x_min(k+1)~=self.x_max(k) && self.t_min(k+1)~=self.t_max(k)))
%                     tempM = zeros(0,0);
%                     %Upper bound constraints
%                     %L1<=M1
%                     %L1<=M2
%                     temp_array = self.mu(k, self.t_min(k+1), self.x_min(k+1));
%                     
%                     for block = 1:size(self.indiced_rho_ini,1) % For initial
%                         
%                         if (self.indiced_rho_ini(block,3) == 0)
%                             
%                             for b = self.indiced_rho_ini(block,1)-1:self.indiced_rho_ini(block,2)-1
%                                 %-1 to be consistent with other functions self.mtau...
%                                 arraym = self.mtau1(b,self.t_min(k+1), self.x_min(k+1));
%                                 array2 = self.substractArray(arraym,temp_array);
%                                 if(~isempty(array2))
%                                     rows = size(list3,1);
%                                     list3(rows+1,:) = array2;
%                                     rows = size(tempM,1);
%                                     tempM(rows+1,:) = arraym;
%                                     break;  % Only need to consider the first one
%                                 end
%                             end
%                             
%                         elseif (self.indiced_rho_ini(block,3) == 1)
%                             
%                             for b = self.indiced_rho_ini(block,2)-1:-1:self.indiced_rho_ini(block,1)-1
%                                 % For k>k_c case, From the last one
%                                 arraym = self.mtau2(b,self.t_min(k+1), self.x_min(k+1));
%                                 array2 = self.substractArray(arraym,temp_array);
%                                 if(~isempty(array2))
%                                     rows = size(list3,1);
%                                     list3(rows+1,:) = array2;
%                                     rows = size(tempM,1);
%                                     tempM(rows+1,:) = arraym;
%                                     break;
%                                 end
%                             end
%                             
%                         elseif (isnan(self.indiced_rho_ini(block,3)))
%                             
%                             for b = self.indiced_rho_ini(block,1)-1:self.indiced_rho_ini(block,2)-1
%                                 % For undefined initial blocks, check all
%                                 arraym = self.mtau1(b,self.t_min(k+1), self.x_min(k+1));
%                                 array2 = self.substractArray(arraym,temp_array);
%                                 if(~isempty(array2))
%                                     rows = size(list3,1);
%                                     list3(rows+1,:) = array2;
%                                     rows = size(tempM,1);
%                                     tempM(rows+1,:) = arraym;
%                                 end
%                                 
%                                 arraym = self.mtau2(b,self.t_min(k+1), self.x_min(k+1));
%                                 array2 = self.substractArray(arraym,temp_array);
%                                 if(~isempty(array2))
%                                     rows = size(list3,1);
%                                     list3(rows+1,:) = array2;
%                                     rows = size(tempM,1);
%                                     tempM(rows+1,:) = arraym;
%                                 end
%                             end
%                             
%                         end
%                     end 
%                     
%                     for n=0:self.n_max    % For Downstream
%                         
%                         arraym = self.mbeta(n,self.t_min(k+1), self.x_min(k+1));
%                         array2 = self.substractArray(arraym,temp_array);
%                         if(~isempty(array2) && arraym(1,self.n_max +1 + n+1) ~= self.Tvec(n+1))
%                             % The second condition ensure we only consider
%                             % the solution at the characteristic domain
%                             %Add a constraint
%                             rows = size(list3,1);
%                             list3(rows+1,:) = array2;
%                             %Assign the solution to temporary Matrix
%                             rows = size(tempM,1);
%                             tempM(rows+1,:) = arraym;
%                             break;  %only one
%                         end
%                         
%                     end
%                     
%                     for n=0:self.n_max    % For upstream
%                         
%                         arraym = self.mgamma(n,self.t_min(k+1), self.x_min(k+1));
%                         array2 = self.substractArray(arraym, temp_array);
%                         if(~isempty(array2) && arraym(1,self.size_row)==0)
%                             %Add a constraint
%                             rows = size(list3,1);
%                             list3(rows+1,:) = array2;
%                             %Assign the solution to temporary Matrix
%                             rows = size(tempM,1);
%                             tempM(rows+1,:) = arraym;
%                             break;
%                         end
%                         
%                     end
%                     
%                     for n=0:self.m_max    % For internal
%                         
%                         if(n~=k)
%                             arraym = self.mmu(n,self.t_min(k+1), self.x_min(k+1));
%                             array2 = self.substractArray(arraym,temp_array);
%                             if(~isempty(array2))
%                                 rows = size(list3,1);
%                                 list3(rows+1,:) = array2;
%                                 rows = size(tempM,1);
%                                 tempM(rows+1,:) = arraym;
%                             end
%                         end
%                         
%                     end
%                     
%                     index_u = find(self.t_u~=self.t0+self.Tcum(self.n_max+2));
%                     for i=1:length(index_u) % For density
%                 
%                         n = index_u(i)-1;
%                         
%                         arraym = self.mups1(n,self.t_min(k+1), self.x_min(k+1));
%                         array2 = self.substractArray(arraym,temp_array);
%                         if(~isempty(array2))
%                             rows = size(list3,1);
%                             list3(rows+1,:) = array2;
%                             rows = size(tempM,1);
%                             tempM(rows+1,:) = arraym;
%                         end
%                         
%                         arraym = self.mups2(n,self.t_min(k+1), self.x_min(k+1));
%                         array2 = self.substractArray(arraym,temp_array);
%                         if(~isempty(array2))
%                             rows = size(list3,1);
%                             list3(rows+1,:) = array2;
%                             rows = size(tempM,1);
%                             tempM(rows+1,:) = arraym;
%                         end
%                         
%                     end
%                     
%                     rowsM = size(tempM,1);
%                     
%                     %Define all the possible combinations according to the solutions that apply
%                     %at the specified point
%                     
%                     %Lower bound constraints
%                     %C*b1 +L1>=M1
%                     %C*b1 +L1>=M2
%                     
%                     for i = 1:rowsM
%                         
%                         array = temp_array;
%                         
%                         %Decode the binary combinations
%                         
%                         for counter=1:self.nb_min(k+1)
%                             if comb_mat(i,nb-self.nb_min(k+1)+counter) == 1
%                                 array(1,2*(self.n_max + 1) + self.b_max + 1 +2*(self.m_max+1) +2*(self.u_max+1)+ countB+ counter) = -Cmax;
%                                 
%                             elseif comb_mat(i,nb-self.nb_min(k+1)+counter) == 0
%                                 array(1,2*(self.n_max + 1) + self.b_max + 1 +2*(self.m_max+1) +2*(self.u_max+1)+ countB+ counter) = Cmax;
%                                 
%                             end
%                         end
%                         array(1,self.size_row) = -Cmax*(sum(comb_mat(i,:))); %RHS (set negative to be on same side)
%                         
%                         % Add constraint to the MILP matrix
%                         
%                         array2 = self.substractArray(array,tempM(i,:));
%                         rows = size(list3,1);
%                         list3(rows+1,:) = array2;
%                     end
%                     
%                     %Define the last constraint to bound the binary combination
%                     
%                     array = zeros(1,self.size_row);
%                     
%                     for counter=1:self.nb_min(k+1)
%                         
%                         array(1,2*(self.n_max + 1) + self.b_max + 1 +2*(self.m_max+1) + 2*(self.u_max+1) + countB + counter) = -2^(self.nb_min(k+1)-counter);
%                         
%                     end
%                     
%                     array(1,self.size_row) = -(rowsM-1); %RHS (maximum possible value of the binary comb)
%                     rows = size(list3,1);
%                     list3(rows+1,:) = array;
%                     
%                     countB = countB+self.nb_min(k+1);
%                     
%                 end
%                 
%                 %x_max and t_max point-------------------------
%                 
%                 tempM = zeros(0,0);
%                 
%                 %Upper bound constraints
%                 %L1<=M1
%                 %L1<=M2
%                 
%                 temp_array = self.mu(k, self.t_max(k+1), self.x_max(k+1));
%                 
%                 for block = 1:size(self.indiced_rho_ini,1) % For initial
%                     
%                     if (self.indiced_rho_ini(block,3) == 0)
%                         
%                         for b = self.indiced_rho_ini(block,1)-1:self.indiced_rho_ini(block,2)-1
%                             %-1 to be consistent with other functions self.mtau...
%                             arraym = self.mtau1(b,self.t_max(k+1), self.x_max(k+1));
%                             array2 = self.substractArray(arraym,temp_array);
%                             if(~isempty(array2))
%                                 rows = size(list3,1);
%                                 list3(rows+1,:) = array2;
%                                 rows = size(tempM,1);
%                                 tempM(rows+1,:) = arraym;
%                                 break;  % Only need to consider the first one
%                             end
%                         end
%                         
%                     elseif (self.indiced_rho_ini(block,3) == 1)
%                         
%                         for b = self.indiced_rho_ini(block,2)-1:-1:self.indiced_rho_ini(block,1)-1
%                             % For k>k_c case, From the last one
%                             arraym = self.mtau2(b,self.t_max(k+1), self.x_max(k+1));
%                             array2 = self.substractArray(arraym,temp_array);
%                             if(~isempty(array2))
%                                 rows = size(list3,1);
%                                 list3(rows+1,:) = array2;
%                                 rows = size(tempM,1);
%                                 tempM(rows+1,:) = arraym;
%                                 break;
%                             end
%                         end
%                         
%                     elseif isnan(self.indiced_rho_ini(block,3))
%                         
%                         for b = self.indiced_rho_ini(block,1)-1:self.indiced_rho_ini(block,2)-1
%                             % For undefined initial blocks, check all
%                             arraym = self.mtau1(b,self.t_max(k+1), self.x_max(k+1));
%                             array2 = self.substractArray(arraym,temp_array);
%                             if(~isempty(array2))
%                                 rows = size(list3,1);
%                                 list3(rows+1,:) = array2;
%                                 rows = size(tempM,1);
%                                 tempM(rows+1,:) = arraym;
%                             end
%                             
%                             arraym = self.mtau2(b,self.t_max(k+1), self.x_max(k+1));
%                             array2 = self.substractArray(arraym,temp_array);
%                             if(~isempty(array2))
%                                 rows = size(list3,1);
%                                 list3(rows+1,:) = array2;
%                                 rows = size(tempM,1);
%                                 tempM(rows+1,:) = arraym;
%                             end
%                         end
%                         
%                     end
%                     
%                 end
%                 
%                 for n=0:self.n_max
%                     arraym = self.mbeta(n,self.t_max(k+1), self.x_max(k+1));
%                     array2 = self.substractArray(arraym,temp_array);
%                     if(~isempty(array2) && array(1,self.n_max+1 + n+1) ~= self.Tvec(n+1))
%                         %Add a constraint
%                         rows = size(list3,1);
%                         list3(rows+1,:) = array2;
%                         %Assign the solution to temporary Matrix
%                         rows = size(tempM,1);
%                         tempM(rows+1,:) = arraym;
%                         break;
%                     end
%                 end
%                 
%                 for n=0:self.n_max
%                     arraym = self.mgamma(n,self.t_max(k+1), self.x_max(k+1));
%                     array2 = self.substractArray(arraym, temp_array);
%                     if(~isempty(array2)  && arraym(1,self.size_row)==0)
%                         %Add a constraint
%                         rows = size(list3,1);
%                         list3(rows+1,:) = array2;
%                         %Assign the solution to temporary Matrix
%                         rows = size(tempM,1);
%                         tempM(rows+1,:) = arraym;
%                         break;
%                     end
%                     
%                 end
%                 
%                 for n=0:self.m_max
%                     
%                     
%                     if(n~=k)
%                         arraym = self.mmu(n,self.t_max(k+1), self.x_max(k+1));
%                         array2 = self.substractArray(arraym,temp_array);
%                         if(~isempty(array2))
%                             %Add a constraint
%                             rows = size(list3,1);
%                             list3(rows+1,:) = array2;
%                             %Assign the solution to temporary Matrix
%                             rows = size(tempM,1);
%                             tempM(rows+1,:) = arraym;
%                             
%                         end
%                         
%                     end
%                     
%                     
%                 end
%                 
%                 index_u = find(self.t_u~=self.t0+self.Tcum(self.n_max+2));
%                 for i=1:length(index_u) % For density
%                     
%                     n = index_u(i)-1;
%                     
%                     arraym = self.mups1(n,self.t_max(k+1), self.x_max(k+1));
%                     array2 = self.substractArray(arraym,temp_array);
%                     if(~isempty(array2))
%                         rows = size(list3,1);
%                         list3(rows+1,:) = array2;
%                         rows = size(tempM,1);
%                         tempM(rows+1,:) = arraym;
%                     end
%                     
%                     arraym = self.mups2(n,self.t_max(k+1), self.x_max(k+1));
%                     array2 = self.substractArray(arraym,temp_array);
%                     if(~isempty(array2))
%                         rows = size(list3,1);
%                         list3(rows+1,:) = array2;
%                         rows = size(tempM,1);
%                         tempM(rows+1,:) = arraym;
%                     end
%                     
%                 end
%                 
%                 rowsM = size(tempM,1);
%                 
%                 %Define all the possible combinations according to the solutions that apply
%                 %at the specified point
%                 
%                 %Lower bound constraints
%                 %C*b1 +L1>=M1
%                 %C*b1 +L1>=M2
%                 
%                 for i = 1:rowsM
%                     
%                     array = temp_array;
%                     
%                     %Decode the binary combinations
%                     for counter=1:self.nb_max(1,k+1)
%                         if comb_mat(i,nb-self.nb_max(1,k+1)+counter) == 1
%                             array(1,2*(self.n_max + 1) + self.b_max + 1 +2*(self.m_max+1) +2*(self.u_max+1)+ countB + counter) = -Cmax;
%                             
%                         elseif comb_mat(i,nb-self.nb_max(1,k+1)+counter) == 0
%                             array(1,2*(self.n_max + 1) + self.b_max + 1 +2*(self.m_max+1) +2*(self.u_max+1)+ countB + counter) = Cmax;
%                             
%                         end
%                     end
%                     array(1,self.size_row) = -Cmax*(sum(comb_mat(i,:))); %RHS (set negative to be on same side)
%                     
%                     % Add constraint to the matrix constraints
%                     
%                     array2 = self.substractArray(array,tempM(i,:));
%                     rows = size(list3,1);
%                     list3(rows+1,:) = array2;
%                 end
%                 
%                 %Define the last constraint to bound the binary combination
%                 
%                 array = zeros(1,self.size_row);
%                 
%                 % Define constraint matrix elements for binary terms
%                 for counter=1:self.nb_max(1,k+1)
%                     
%                     array(1,2*(self.n_max + 1) + self.b_max + 1 +2*(self.m_max+1) + 2*(self.u_max+1)+ countB + counter) = -2^(self.nb_max(1,k+1)-counter);
%                     
%                 end
%                 
%                 array(1,self.size_row) = -(rowsM-1); %RHS
%                 rows = size(list3,1);
%                 list3(rows+1,:) = array;
%                 
%                 countB = countB + self.nb_max(k+1);
%             end
% 
%         end
%         
%         
%         function [list3] = setMatrix3u(si)
%             
%             list3 = zeros(0,0);
%             
%             %Binary variables useful counter as a reference
%             countB=0;
%             Cmax = 500000; 
%             
%             % nb: number of binary variables per point
%             nb = ceil(log2(2*(self.b_max+1)+2*(self.n_max+1)+(self.m_max+1)+(self.u_max+1)));
%             
%             % comb_mat: combination matrix for possible binary combinations
%             comb_mat = zeros(2^nb,nb);
%             
%             for i = 1:2^nb
%                 comb_mat(i,:) = double(dec2bin(i-1,nb)=='1');
%             end
%             
%             % Binary inequalities induced by density conditions
%             for k = 0:self.u_max
%                  %x_min_u and t_u point------------------
%                 if (k==0 || (self.x_min_u(k+1)~=self.x_max_u(k) && self.t_u(k+1)~=self.t_u(k)))
%                     tempM = zeros(0,0);
%                     %Upper bound constraints
%                     %L1<=M1
%                     %L1<=M2
%                     temp_array = self.ups(k, self.t_u(k+1), self.x_min_u(k+1));
%                     
%                     for block = 1:size(self.indiced_rho_ini,1) % For initial
%                         
%                         if (self.indiced_rho_ini(block,3) == 0)
%                             
%                             for b = self.indiced_rho_ini(block,1)-1:self.indiced_rho_ini(block,2)-1
%                                 %-1 to be consistent with other functions self.mtau...
%                                 arraym = self.mtau1(b,self.t_u(k+1), self.x_min_u(k+1));
%                                 array2 = self.substractArray(arraym,temp_array);
%                                 if(~isempty(array2))
%                                     rows = size(list3,1);
%                                     list3(rows+1,:) = array2;
%                                     rows = size(tempM,1);
%                                     tempM(rows+1,:) = arraym;
%                                     break;  % Only need to consider the first one
%                                 end
%                             end
%                             
%                         elseif (self.indiced_rho_ini(block,3) == 1)
%                             
%                             for b = self.indiced_rho_ini(block,2)-1:-1:self.indiced_rho_ini(block,1)-1
%                                 % For k>k_c case, From the last one
%                                 arraym = self.mtau2(b,self.t_u(k+1), self.x_min_u(k+1));
%                                 array2 = self.substractArray(arraym,temp_array);
%                                 if(~isempty(array2))
%                                     rows = size(list3,1);
%                                     list3(rows+1,:) = array2;
%                                     rows = size(tempM,1);
%                                     tempM(rows+1,:) = arraym;
%                                     break;
%                                 end
%                             end
%                             
%                         elseif (isnan(self.indiced_rho_ini(block,3)))
%                             
%                             for b = self.indiced_rho_ini(block,1)-1:self.indiced_rho_ini(block,2)-1
%                                 % For undefined initial blocks, check all
%                                 arraym = self.mtau1(b,self.t_u(k+1), self.x_min_u(k+1));
%                                 array2 = self.substractArray(arraym,temp_array);
%                                 if(~isempty(array2))
%                                     rows = size(list3,1);
%                                     list3(rows+1,:) = array2;
%                                     rows = size(tempM,1);
%                                     tempM(rows+1,:) = arraym;
%                                 end
%                                 
%                                 arraym = self.mtau2(b,self.t_u(k+1), self.x_min_u(k+1));
%                                 array2 = self.substractArray(arraym,temp_array);
%                                 if(~isempty(array2))
%                                     rows = size(list3,1);
%                                     list3(rows+1,:) = array2;
%                                     rows = size(tempM,1);
%                                     tempM(rows+1,:) = arraym;
%                                 end
%                             end
%                             
%                         end
%                     end 
%                     
%                     for n=0:self.n_max    % For Downstream
%                         
%                         arraym = self.mbeta(n,self.t_u(k+1), self.x_min_u(k+1));
%                         array2 = self.substractArray(arraym,temp_array);
%                         if(~isempty(array2) && arraym(1,self.n_max + n+2) ~= self.Tvec(n+1))
%                             %Add a constraint
%                             rows = size(list3,1);
%                             list3(rows+1,:) = array2;
%                             %Assign the solution to temporary Matrix
%                             rows = size(tempM,1);
%                             tempM(rows+1,:) = arraym;
%                             break;  %only one
%                         end
%                         
%                     end
%                     
%                     for n=0:self.n_max    % For upstream
%                         
%                         arraym = self.mgamma(n,self.t_u(k+1), self.x_min_u(k+1));
%                         array2 = self.substractArray(arraym, temp_array);
%                         if(~isempty(array2) && arraym(1,self.size_row)==0)
%                             %Add a constraint
%                             rows = size(list3,1);
%                             list3(rows+1,:) = array2;
%                             %Assign the solution to temporary Matrix
%                             rows = size(tempM,1);
%                             tempM(rows+1,:) = arraym;
%                             break;
%                         end
%                         
%                     end
%                     
%                     for n=0:self.m_max    % For internal
%                         
%                         arraym = self.mmu(n,self.t_u(k+1), self.x_min_u(k+1));
%                         array2 = self.substractArray(arraym,temp_array);
%                         if(~isempty(array2))
%                             rows = size(list3,1);
%                             list3(rows+1,:) = array2;
%                             rows = size(tempM,1);
%                             tempM(rows+1,:) = arraym;
%                         end
%                         
%                     end
%                     
%                     index_u = find(self.t_u <= self.t0+self.Tcum(self.n_max+2));
%                     for i=1:length(index_u) % For density
%                 
%                         n = index_u(i)-1;
%                         
%                         
%                             arraym = self.mups1(n,self.t_u(k+1), self.x_min_u(k+1));
%                             array2 = self.substractArray(arraym,temp_array);
%                             if(~isempty(array2))
%                                 rows = size(list3,1);
%                                 list3(rows+1,:) = array2;
%                                 rows = size(tempM,1);
%                                 tempM(rows+1,:) = arraym;
%                             end
%                             
%                             arraym = self.mups2(n,self.t_u(k+1), self.x_min_u(k+1));
%                             array2 = self.substractArray(arraym,temp_array);
%                             if(~isempty(array2))
%                                 rows = size(list3,1);
%                                 list3(rows+1,:) = array2;
%                                 rows = size(tempM,1);
%                                 tempM(rows+1,:) = arraym;
%                             end
%                        
%                         
%                     end
%                     
%                     rowsM = size(tempM,1);
%                     
%                     %Define all the possible combinations according to the solutions that apply
%                     %at the specified point
%                     
%                     %Lower bound constraints
%                     %C*b1 +L1>=M1
%                     %C*b1 +L1>=M2
%                     
%                     for i = 1:rowsM
%                         
%                         array = temp_array;
%                         
%                         %Decode the binary combinations
%                         
%                         for counter=1:self.nb_min_u(k+1)
%                             if comb_mat(i,nb-self.nb_min_u(k+1)+counter) == 1
%                                 array(1,2*(self.n_max + 1) + self.b_max + 1 +2*(self.m_max+1) +2*(self.u_max+1) +sum(self.nb_min) +sum(self.nb_max) + countB+ counter) = -Cmax;
%                                 
%                             elseif comb_mat(i,nb-self.nb_min_u(k+1)+counter) == 0
%                                 array(1,2*(self.n_max + 1) + self.b_max + 1 +2*(self.m_max+1) +2*(self.u_max+1)+sum(self.nb_min) +sum(self.nb_max)+ countB+ counter) = Cmax;
%                                 
%                             end
%                         end
%                         array(1,self.size_row) = -Cmax*(sum(comb_mat(i,:))); %RHS (set negative to be on same side)
%                         
%                         % Add constraint to the MILP matrix
%                         
%                         array2 = self.substractArray(array,tempM(i,:));
%                         rows = size(list3,1);
%                         list3(rows+1,:) = array2;
%                     end
%                     
%                     %Define the last constraint to bound the binary combination
%                     array = zeros(1,self.size_row);
%                     
%                     for counter=1:self.nb_min_u(k+1)
%                         
%                         array(1,2*(self.n_max + 1) + self.b_max + 1 +2*(self.m_max+1) +2*(self.u_max+1)+sum(self.nb_min)+sum(self.nb_max)+ countB + counter) = -2^(self.nb_min_u(k+1)-counter);
%                         
%                     end
%                     
%                     array(1,self.size_row) = -(rowsM-1); %RHS (maximum possible value of the binary comb)
%                     rows = size(list3,1);
%                     list3(rows+1,:) = array;
%                     
%                     countB = countB+self.nb_min_u(k+1);
%                     
%                 end
%                 
%                 %x_max_u and t_u point-------------------------
%                 
%                 tempM = zeros(0,0);
%                 
%                 %Upper bound constraints
%                 %L1<=M1
%                 %L1<=M2
%                 
%                 temp_array = self.ups(k, self.t_u(k+1), self.x_max_u(k+1));
%                 
%                 for block = 1:size(self.indiced_rho_ini,1) % For initial
%                     
%                     if (self.indiced_rho_ini(block,3) == 0)
%                         
%                         for b = self.indiced_rho_ini(block,1)-1:self.indiced_rho_ini(block,2)-1
%                             %-1 to be consistent with other functions self.mtau...
%                             arraym = self.mtau1(b,self.t_u(k+1), self.x_max_u(k+1));
%                             array2 = self.substractArray(arraym,temp_array);
%                             if(~isempty(array2))
%                                 rows = size(list3,1);
%                                 list3(rows+1,:) = array2;
%                                 rows = size(tempM,1);
%                                 tempM(rows+1,:) = arraym;
%                                 break;  % Only need to consider the first one
%                             end
%                         end
%                         
%                     elseif (self.indiced_rho_ini(block,3) == 1)
%                         
%                         for b = self.indiced_rho_ini(block,2)-1:-1:self.indiced_rho_ini(block,1)-1
%                             % For k>k_c case, From the last one
%                             arraym = self.mtau2(b,self.t_u(k+1), self.x_max_u(k+1));
%                             array2 = self.substractArray(arraym,temp_array);
%                             if(~isempty(array2))
%                                 rows = size(list3,1);
%                                 list3(rows+1,:) = array2;
%                                 rows = size(tempM,1);
%                                 tempM(rows+1,:) = arraym;
%                                 break;
%                             end
%                         end
%                         
%                     elseif isnan(self.indiced_rho_ini(block,3))
%                         
%                         for b = self.indiced_rho_ini(block,1)-1:self.indiced_rho_ini(block,2)-1
%                             % For undefined initial blocks, check all
%                             arraym = self.mtau1(b,self.t_u(k+1), self.x_max_u(k+1));
%                             array2 = self.substractArray(arraym,temp_array);
%                             if(~isempty(array2))
%                                 rows = size(list3,1);
%                                 list3(rows+1,:) = array2;
%                                 rows = size(tempM,1);
%                                 tempM(rows+1,:) = arraym;
%                             end
%                             
%                             arraym = self.mtau2(b,self.t_u(k+1), self.x_max_u(k+1));
%                             array2 = self.substractArray(arraym,temp_array);
%                             if(~isempty(array2))
%                                 rows = size(list3,1);
%                                 list3(rows+1,:) = array2;
%                                 rows = size(tempM,1);
%                                 tempM(rows+1,:) = arraym;
%                             end
%                         end
%                         
%                     end
%                     
%                 end
%                 
%                 for n=0:self.n_max
%                     arraym = self.mbeta(n,self.t_u(k+1), self.x_max_u(k+1));
%                     array2 = self.substractArray(arraym,temp_array);
%                     if(~isempty(array2) && array(1,self.n_max + n+2) ~= self.Tvec(n+1))
%                         %Add a constraint
%                         rows = size(list3,1);
%                         list3(rows+1,:) = array2;
%                         %Assign the solution to temporary Matrix
%                         rows = size(tempM,1);
%                         tempM(rows+1,:) = arraym;
%                         break;
%                     end
%                 end
%                 
%                 for n=0:self.n_max
%                     arraym = self.mgamma(n,self.t_u(k+1), self.x_max_u(k+1));
%                     array2 = self.substractArray(arraym, temp_array);
%                     if(~isempty(array2)  && arraym(1,self.size_row)==0)
%                         %Add a constraint
%                         rows = size(list3,1);
%                         list3(rows+1,:) = array2;
%                         %Assign the solution to temporary Matrix
%                         rows = size(tempM,1);
%                         tempM(rows+1,:) = arraym;
%                         break;
%                     end
%                     
%                 end
%                 
%                 for n=0:self.m_max
%                     
%                     arraym = self.mmu(n,self.t_u(k+1), self.x_max_u(k+1));
%                     array2 = self.substractArray(arraym,temp_array);
%                     if(~isempty(array2))
%                         %Add a constraint
%                         rows = size(list3,1);
%                         list3(rows+1,:) = array2;
%                         %Assign the solution to temporary Matrix
%                         rows = size(tempM,1);
%                         tempM(rows+1,:) = arraym;
%                     end
%                     
%                 end
%                 
%                 index_u = find(self.t_u <= self.t0+self.Tcum(self.n_max+2));
%                 for i=1:length(index_u) % For density
%                     
%                     n = index_u(i)-1;
%                     
%                
%                         arraym = self.mups1(n,self.t_u(k+1), self.x_max_u(k+1));
%                         array2 = self.substractArray(arraym,temp_array);
%                         if(~isempty(array2))
%                             rows = size(list3,1);
%                             list3(rows+1,:) = array2;
%                             rows = size(tempM,1);
%                             tempM(rows+1,:) = arraym;
%                         end
%                         
%                         arraym = self.mups2(n,self.t_u(k+1), self.x_max_u(k+1));
%                         array2 = self.substractArray(arraym,temp_array);
%                         if(~isempty(array2))
%                             rows = size(list3,1);
%                             list3(rows+1,:) = array2;
%                             rows = size(tempM,1);
%                             tempM(rows+1,:) = arraym;
%                         end
%                  
%                     
%                 end
%                 
%                 rowsM = size(tempM,1);
%                 
%                 %Define all the possible combinations according to the solutions that apply
%                 %at the specified point
%                 
%                 %Lower bound constraints
%                 %C*b1 +L1>=M1
%                 %C*b1 +L1>=M2
%                 
%                 for i = 1:rowsM
%                     
%                     array = temp_array;
%                     
%                     %Decode the binary combinations
%                     for counter=1:self.nb_max_u(1,k+1)
%                         if comb_mat(i,nb-self.nb_max_u(1,k+1)+counter) == 1
%                             array(1,2*(self.n_max + 1) + self.b_max + 1 +2*(self.m_max+1) +  2*(self.u_max+1) + sum(self.nb_min) + sum(self.nb_max) +countB + counter) = -Cmax;
%                             
%                         elseif comb_mat(i,nb-self.nb_max_u(1,k+1)+counter) == 0
%                             array(1,2*(self.n_max + 1) + self.b_max + 1 +2*(self.m_max+1) + 2*(self.u_max+1) + sum(self.nb_min) + sum(self.nb_max) + countB + counter) = Cmax;
%                             
%                         end
%                     end
%                     array(1,self.size_row) = -Cmax*(sum(comb_mat(i,:))); %RHS (set negative to be on same side)
%                     
%                     % Add constraint to the matrix constraints
%                     
%                     array2 = self.substractArray(array,tempM(i,:));
%                     rows = size(list3,1);
%                     list3(rows+1,:) = array2;
%                 end
%                 
%                 %Define the last constraint to bound the binary combination
%                 
%                 array = zeros(1,self.size_row);
%                 
%                 % Define constraint matrix elements for binary terms
%                 for counter=1:self.nb_max_u(1,k+1)
%                     
%                     array(1,2*(self.n_max + 1) + self.b_max + 1 +2*(self.m_max+1)+ 2*(self.u_max+1) + sum(self.nb_min) + sum(self.nb_max) + countB + counter) = -2^(self.nb_max_u(1,k+1)-counter);
%                     
%                 end
%                 
%                 array(1,self.size_row) = -(rowsM-1); %RHS
%                 rows = size(list3,1);
%                 list3(rows+1,:) = array;
%                 
%                 countB = countB + self.nb_max_u(k+1);
%           
%             end
%             
%         end
%         
%         
%         function [list4] = setMatrixAux(si)
%              list4 = zeros(0,0);
%              
%              for u=0:self.u_max
%                  array = zeros(1,self.size_row);
%                  array(1,2*(self.n_max+1) + (self.b_max+1) + 2*(self.m_max + 1) + (self.u_max+1) + u + 1) = -1;
%                  array(1,2*(self.n_max+1) + (self.b_max+1) + 2*(self.m_max + 1) + 2*(self.u_max+1) + sum(self.nb_min) + sum(self.nb_max)...
%                      + sum(self.nb_min_u) + sum(self.nb_max_u) + u + 1) = self.k_m-self.k_c;
%                  array(1,self.size_row) = -self.k_c;
%                  rows = size(list4,1);
%                  list4(rows+1,:) = array;
%              end
%              
%              for u=0:self.u_max
%                  array = zeros(1,self.size_row);
%                  array(1,2*(self.n_max+1) + (self.b_max+1) + 2*(self.m_max + 1) + (self.u_max+1) + u + 1) = 1;
%                  array(1,2*(self.n_max+1) + (self.b_max+1) + 2*(self.m_max + 1) + 2*(self.u_max+1) + sum(self.nb_min) + sum(self.nb_max)...
%                      + sum(self.nb_min_u) + sum(self.nb_max_u) + u + 1) = -self.k_c;
%                  array(1,self.size_row) = 0;
%                  rows = size(list4,1);
%                  list4(rows+1,:) = array;
%              end
%             
%         end
%             
        
        
        %=========================================================================
        %Get Number Binary variables
        
        function [nb_min,nb_max, nb_min_u, nb_max_u] = getBinaryvar(self)
            
            nb_min = zeros(1,self.m_max);
            nb_max = zeros(1,self.m_max);
            
            nb_min_u = zeros(1,self.u_max);
            nb_max_u = zeros(1,self.u_max);
            
            % Binaries induced by internal conditions
            for k = 1:self.m_max
                % For x_min
                if (k==1 || (self.x_min(k)~=self.x_max(k-1) &&...
                        self.t_min(k)~=self.t_max(k-1)))
                    
                    countM=0;
                    
                    % For initial
                    for block = 1:size(self.indiced_rho_ini,1) 
                        
                        if (self.indiced_rho_ini(block,3) == 0)   % For k<k_c case, only the first one
                            
                            for b = self.indiced_rho_ini(block,1):self.indiced_rho_ini(block,2)
                                array = self.mtau1(b,self.t_min(k), self.x_min(k));
                                if(~isempty(array))
                                    countM = countM+1;
                                    break;  
                                end
                            end
                            
                        elseif (self.indiced_rho_ini(block,3) == 1) % For k>k_c case, only the last one
                            
                            for b = self.indiced_rho_ini(block,2):-1:self.indiced_rho_ini(block,1)
                                array = self.mtau2(b,self.t_min(k), self.x_min(k));
                                if(~isempty(array))
                                    countM = countM+1;
                                    break;
                                end
                            end
                            
                        elseif isnan(self.indiced_rho_ini(block,3))
                            
                            for b = self.indiced_rho_ini(block,1):self.indiced_rho_ini(block,2)
                                array = self.mtau1(b,self.t_min(k), self.x_min(k));
                                if(~isempty(array))
                                    countM = countM+1;
                                end
                                
                                array = self.mtau2(b,self.t_min(k), self.x_min(k));
                                if(~isempty(array))
                                    countM = countM+1;
                                end
                            end
                            
                        end
                        
                    end
                    
                    % for downstream
                    % in downstream characteristic domain
                    for n=1:self.n_max_ds
                        
                        array = self.mbeta(n,self.t_min(k), self.x_min(k));
                        if(~isempty(array) && array(1,self.n_max_us + n) ~= self.T_ds(n))
                            % If in characteristic domain
                            countM = countM+1;
                            break;
                        end
                        
                    end
                    
                    % for upstream
                    % in upstream characteristic domain
                    for n=1:self.n_max_us
                        
                        array = self.mgamma(n,self.t_min(k), self.x_min(k));
                        if(~isempty(array) && array(1,n) ~= self.T_us(n) )
                            countM = countM+1;
                            break;
                        end
                        
                    end
                    
                    % for other internal conditions
                    for n=1:self.m_max
                        
                        if (n~=k)
                            array = self.mmu(n,self.t_min(k), self.x_min(k));
                            if(~isempty(array))
                                countM = countM+1;
                            end
                        end
                        
                    end
                    
                    % for density conditions
                    index_u = find(self.t_u~= self.end_time);
                    for i=1:length(index_u) % For density
                
                        n = index_u(i);
                        
                        array = self.mups1(n,self.t_min(k), self.x_min(k));
                        if(~isempty(array))
                            countM = countM+1;
                        end
                        
                        array = self.mups2(n,self.t_min(k), self.x_min(k));
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
                            array = self.mtau1(b,self.t_min(k), self.x_min(k));
                            if(~isempty(array))
                                countM = countM+1;
                                break;
                            end
                        end
                        
                    elseif (self.indiced_rho_ini(block,3) == 1) % For k>k_c case, only the last one
                        
                        for b = self.indiced_rho_ini(block,2):-1:self.indiced_rho_ini(block,1)
                            array = self.mtau2(b,self.t_min(k), self.x_min(k));
                            if(~isempty(array))
                                countM = countM+1;
                                break;
                            end
                        end
                        
                    elseif isnan(self.indiced_rho_ini(block,3))
                        
                        for b = self.indiced_rho_ini(block,1):self.indiced_rho_ini(block,2)
                            array = self.mtau1(b,self.t_min(k), self.x_min(k));
                            if(~isempty(array))
                                countM = countM+1;
                            end
                            
                            array = self.mtau2(b,self.t_min(k), self.x_min(k));
                            if(~isempty(array))
                                countM = countM+1;
                            end
                        end
                        
                    end
                    
                end
                
                % for downstream
                % in downstream characteristic domain
                for n=1:self.n_max_ds
                    
                    array = self.mbeta(n,self.t_min(k), self.x_min(k));
                    if(~isempty(array) && array(1,self.n_max_us + n) ~= self.T_ds(n))
                        % If in characteristic domain
                        countM = countM+1;
                        break;
                    end
                    
                end
                
                % for upstream
                % in upstream characteristic domain
                for n=1:self.n_max_us
                    
                    array = self.mgamma(n,self.t_min(k), self.x_min(k));
                    if(~isempty(array) && array(1,n) ~= self.T_us(n) )
                        countM = countM+1;
                        break;
                    end
                    
                end
                
                % for other internal conditions
                for n=0:self.m_max
                    
                    if (n~=k)
                        array = self.mmu(n,self.t_min(k), self.x_min(k));
                        if(~isempty(array))
                            countM = countM+1;
                        end
                    end
                    
                end
                
                % for density conditions
                index_u = find(self.t_u~= self.end_time);
                for i=1:length(index_u) % For density
                    
                    n = index_u(i);
                    
                    array = self.mups1(n,self.t_min(k), self.x_min(k));
                    if(~isempty(array))
                        countM = countM+1;
                    end
                    
                    array = self.mups2(n,self.t_min(k), self.x_min(k));
                    if(~isempty(array))
                        countM = countM+1;
                    end
                    
                end
                
                nb_max(1,k) = ceil(log2(countM));
            end
            
            
            %Binaries induced by density conditions
            for k = 1:self.u_max
                % For x_min
                if (k==1 || (self.x_min_u(k)~=self.x_max_u(k-1) && self.t_u(k)~=self.t_u(k-1)))
                    
                    countM=0;
                    
                    for block = 1:size(self.indiced_rho_ini,1) % For initial
                        
                        if (self.indiced_rho_ini(block,3) == 0)   % For k<k_c case, only the first one
                            
                            for b = self.indiced_rho_ini(block,1):self.indiced_rho_ini(block,2)
                                array = self.mtau1(b,self.t_u(k), self.x_min_u(k));
                                if(~isempty(array))
                                    countM = countM+1;
                                    break;  
                                end
                            end
                            
                        elseif (self.indiced_rho_ini(block,3) == 1) % For k>k_c case, only the last one
                            
                            for b = self.indiced_rho_ini(block,2):-1:self.indiced_rho_ini(block,1)
                                array = self.mtau2(b,self.t_u(k), self.x_min_u(k));
                                if(~isempty(array))
                                    countM = countM+1;
                                    break;
                                end
                            end
                            
                        elseif isnan(self.indiced_rho_ini(block,3))
                            
                            for b = self.indiced_rho_ini(block,1):self.indiced_rho_ini(block,2)
                                array = self.mtau1(b,self.t_u(k), self.x_min_u(k));
                                if(~isempty(array))
                                    countM = countM+1;
                                end
                                
                                array = self.mtau2(b,self.t_u(k), self.x_min_u(k));
                                if(~isempty(array))
                                    countM = countM+1;
                                end
                            end
                            
                        end
                        
                    end
                    
                    % downstream
                    for n=1:self.n_max_ds
                        
                        array = self.mbeta(n,self.t_u(k), self.x_min_u(k));
                        if(~isempty(array) && array(1,self.n_max_us + n) ~= self.T_ds(n))
                            countM = countM+1;
                            break;
                        end
                        
                    end
                    
                    % upstream
                    for n=1:self.n_max_us
                        
                        array = self.mgamma(n,self.t_u(k), self.x_min_u(k));
                        if(~isempty(array) && array(1,n) ~= self.T_us(n) )
                            countM = countM+1;
                            break;
                        end
                        
                    end
                    
                    % internal conditions
                    for n=1:self.m_max
                        
                        array = self.mmu(n,self.t_u(k), self.x_min_u(k));
                        if(~isempty(array))
                            countM = countM+1;
                        end
                        
                    end
                    
                    index_u = find(self.t_u~= self.end_time);
                    for i=1:length(index_u) % For density
                
                        n = index_u(i);
                        
                  
                            array = self.mups1(n,self.t_u(k), self.x_min_u(k));
                            if(~isempty(array))
                                countM = countM+1;
                            end
                            
                            array = self.mups2(n,self.t_u(k), self.x_min_u(k));
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
                            array = self.mtau1(b,self.t_u(k), self.x_max_u(k));
                            if(~isempty(array))
                                countM = countM+1;
                                break;  % Only need to consider the first one
                            end
                        end
                        
                    elseif (self.indiced_rho_ini(block,3) == 1)
                        
                        for b = self.indiced_rho_ini(block,2):-1:self.indiced_rho_ini(block,1)
                            % For k>k_c case, From the last one
                            array = self.mtau2(b,self.t_u(k), self.x_max_u(k));
                            if(~isempty(array))
                                countM = countM+1;
                                break;
                            end
                        end
                        
                    elseif isnan(self.indiced_rho_ini(block,3))
                        
                        for b = self.indiced_rho_ini(block,1)-1:self.indiced_rho_ini(block,2)-1
                            array = self.mtau1(b,self.t_u(k), self.x_max_u(k));
                            if(~isempty(array))
                                countM = countM+1;
                            end
                            
                            array = self.mtau2(b,self.t_u(k), self.x_max_u(k));
                            if(~isempty(array))
                                countM = countM+1;
                            end
                        end
                        
                    end
                end
                
                % downstream
                for n=1:self.n_max_ds
                    array = self.mbeta(n,self.t_u(k), self.x_max_u(k));
                    if(~isempty(array) && array(1,self.n_max_us + n) ~= self.T_ds(n))
                        countM = countM+1;
                        break;
                    end
                end
                
                % upstream
                for n=1:self.n_max_us
                    array = self.mgamma(n,self.t_u(k), self.x_max_u(k));
                    if(~isempty(array) && array(1,n) ~= self.T_us(n) )
                        countM = countM+1;
                        break;
                    end
                end
                
                % internal
                for n=1:self.m_max
                    
                        array = self.mmu(n,self.t_u(k), self.x_max_u(k));
                        if(~isempty(array))
                            countM = countM+1;
                        end
              
                end
                
                index_u = find(self.t_u ~= self.end_time);
                for i=1:length(index_u) % For density
                    
                    n = index_u(i);
                    
                 
                        array = self.mups1(n,self.t_u(k), self.x_max_u(k));
                        if(~isempty(array))
                            countM = countM+1;
                        end
                        
                        array = self.mups2(n,self.t_u(k), self.x_max_u(k));
                        if(~isempty(array))
                            countM = countM+1;
                        end
               
                    
                end
                
                nb_max_u(1,k) = ceil(log2(countM));
            end
      
        end
        
        
        %===============================================================
        % utility function
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
    
        
        
        
    end
end





