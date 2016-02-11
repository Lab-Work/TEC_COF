% Yanning Li, sep 07, 2015
% This class is for set up a model predictive controller
% It can be used for onramp metering of offramp control. For now, we only
% consider the work zone control scenario with only one merge or one
% diverge.
% The control scheme is as follows::
% 1. Model predictive control scheme. E.g. use past period boundary measurement 
%    data and the historical data for the next period, then output optimal
%    control in the predict period.
% 2. Roll forward once new measurement is available. The rolling step can
%    be 5 min.
% 3. Only the first period of the 30 min optimal control output will be
%    applied on the ramp. Once new optimal control is availabel (updated by
%    new measurement), the new optimal will be applied to the ramp.
% 4. This class simulates a real-time control. The computation time of the
%    solver is taken into account. Only the optimal control output after
%    t_now + dt_computation time will be applied.
% REMARK: 
% - All communication between MATLAB and AIMSUN are in veh/hr, veh/km,
%   km/hr, and s
% - MATLAB internally use m and s. 
% - AIMSUN internally may use km/hr or mile/hr depending on configuration.



classdef rampController < handle
    
    properties
        
        % The network to be controlled
        net;
        
        % historical data
        historical_data;
        
        % past period data
        past_period_data;
        
        % predict period data (extracted from historical data)
        predict_period_data;
        
        % all measurement data
        all_measurement_data;
        
        % all discretized control signal.
        % [setting_time, duration, flow (veh/s)]
        all_signal;
        
        % the short period signal to be write to file and applied in AIMSUN
        % [setting_time, duration, flow (veh/s)]
        % When writing to file, veh/s was converted to veh/hr
        signal_to_apply;
        
        % the starting time of the entire simulation. This is the abusolute
        % zero time of the system.
        t_horizon_start;
        
        % the end time of the entire simulation. The simulation will stop
        % at this time.
        t_horizon_end;
        
%         % rolling step duration. This parameter should be the interval
%         % between two measurement. Set this parameter to prevent updating
%         % the control signal too fast
%         dt_roll_step;
        
        % the duration for warming up.
        dt_warm_up;
        
        % Use the past dt_past measurement data to predict the future.
        % This value must be larger than dt_warm_up
        dt_past;
        
        % compute the optimal control for dt_predict in the future
        dt_predict;
        
        % the current time in system
        t_now;
        
        % the following two are the simulated time for CP
        % t_sim_start = max(t_horizon_start, t_now-dt_past)
        % t_sim_end = min(t_horizon_edn, t_now+dt_predict)
        t_sim_start;
        t_sim_end;
        
        % other paramters that can be leave as default
        % the entropic solution tolerance in veh
        entropy_tolerance;  
        
        % The time space diagram resolution for visualization
        dx_res;     % m
        dt_res;     % s
        
        % the type of the controller: 'none', 'onramp', 'offramp'
        control_type;
        
%         % accuracy of different types of data
%         % The accuracy of the historical data
%         e_historical_flow; 
%         
%         % The accuracy of new boundary flow data
%         e_flow_measurement;
%         e_speed_measurement;
%         
%         % The estimation error of the initial density
%         e_init_estimation;
        
        % The struct that saves the communication file paths
        com;
        % it contains:
        % - the network file written by AIMSUN, containing the network
        %   topology and properties
        %   net_write_done;
        %   net;
        % - the new boundary data file written by AIMSUN
        %   data_write_done;
        %   data;
        % - the updated signal file written by MATLAB
        %   signal_write_done;  % if exists, meaning signal updated
        %   signal;
        % - the initialization of MATLAB is done
        %   matlab_init_done;
        % - the flag for the simulation completed
        %   simulation_completed;
        
    end
    
    
    
    methods
        
        %===============================================================
        % initialize the controller
        function self = rampController(~)
            
            self.t_horizon_start = 0;
            self.t_now = 0;
            
            self.historical_data = struct;
            self.past_period_data = struct;
            self.predict_period_data = struct;
            
            % some that can be leave as default
            self.entropy_tolerance = 1;  
            
            self.dx_res = 5;    % m
            self.dt_res = 2;    % s
            
%             self.e_historical_flow = 0.3;
%             self.e_flow_measurement = 0.05;
%             self.e_speed_measurement = 0.1;
%             
%             self.e_init_estimation = 0.1;
            
            self.net = struct;
            
        end
        
        
        %===============================================================
        % configure the control parameters
        function configController(self, simulation_start_time,... 
                                  simulation_end_time,...
                                  dt_past, dt_predict, ...
                                  control_type)
                              
            self.t_horizon_start = simulation_start_time;               
            self.t_horizon_end = simulation_end_time;
            self.dt_past = dt_past;
            self.dt_predict = dt_predict;
%             self.dt_roll_step = dt_roll_step;
            
            self.control_type = control_type;
            
        end
        
        
        %===============================================================
        % set up the communication between MATLAB controller and AIMSUN
        % simulator using file shareing
        function setUpCommunication(self, com_config_file)
            
            % MATLAB and AIMSUN communicate using shared files
            % set up the file path and names in the COM_CONFIG.m file
            com_config = self.readTextFile(com_config_file);
            
            for line = com_config
                
                % parse and set corresponding paths
                items = strsplit(line{1},',');
                
                self.com.(items{1}) = items{2};
%                 if strcmp(items{1}, 'net')
%                     self.com_net_file = items{2};
%                 elseif strcmp(items{1}, 'net_write_done')
%                     self.com_net_write_done = items{2};
%                 elseif strcmp(items{1}, 'data')
%                     self.com_data_file = items{2};
%                 elseif strcmp(items{1}, 'data_write_done')
%                     self.com_data_write_done = items{2};
%                 elseif strcmp(items{1}, 'signal')
%                     self.com_signal_file = items{2};
%                 elseif strcmp(items{1}, 'signal_write_done')
%                     self.com_signal_write_done = items{2};
%                 elseif strcmp(items{1}, 'matlab_init_done')
%                     self.com_matlab_init_done = items{2};
%                 end
            end
            
        end
        
        
        
        %===============================================================
        % read network configuration
        % 1. read the net work configuration from file, auto generated by 
        %    AIMSUN. 
        % 2. Create a network structure net in this object
        function readNetworkFile(self)
            
            % wait for the network parameter files from CORSIM
            disp('Waiting for AIMSUN network information...');
            while (~exist(self.com.net_write_done,'file'))
                pause(0.1);     % wait 0.5 second
            end
            disp('Got AIMSUN network information...');
            
            % set up the network
            self.net = initNetwork();
            
            % parse and set up the network.
            content = self.readTextFile(self.com.net);
            for line = content
                
                items = strsplit(line{1}, ';');
                
                % if it is a link
                if strcmp(items{1},'link')
                    
                    % build a parameter struct from file
                    % all unit in meters and seconds, exept lengthKM in KM
                    tmp_para = struct;
                    
                    % parse each property
                    for i = 2:length(items)
                        
                        tup = strsplit(items{i}, ',');
                        
                        if strcmp(tup{1},'link_id')
                            link_id = str2double(tup{2});
                        elseif strcmp(tup{1}, 'num_lanes')
                            num_lanes = str2double(tup{2});
                        elseif strcmp(tup{1}, 'lengthKM')
                            len_km = str2double(tup{2});
                        elseif strcmp(tup{1}, 'link_type')
                            link_type = tup{2};
                        elseif strcmp(tup{1}, 'vf')
                            tmp_para.vf = str2double(tup{2});
                        elseif strcmp(tup{1}, 'w')
                            tmp_para.w = str2double(tup{2});   
                        elseif strcmp(tup{1}, 'kc_pl')
                            tmp_para.kc_pl = str2double(tup{2});    
                        elseif strcmp(tup{1}, 'km_pl')
                            tmp_para.km_pl = str2double(tup{2});   
                        elseif strcmp(tup{1}, 'qmax_pl')
                            tmp_para.qmax_pl = str2double(tup{2});
                        elseif strcmp(tup{1}, 'v_min')
                            tmp_para.v_min = str2double(tup{2});
                        elseif strcmp(tup{1}, 'v_max')
                            tmp_para.v_max = str2double(tup{2});
                        end
                    end
                    
                     % initialize the link in the net
                    self.net.addLink( link_id,tmp_para,...
                                          num_lanes,len_km,...
                                          link_type)
                    
                elseif strcmp(items{1},'junc')    
                    % if it is a junction
                    % parse each property
                    for i = 2:length(items)
                        
                        tup = strsplit(items{i}, ',');
                        
                        if strcmp(tup{1}, 'junc_id')
                            junc_id = str2double(tup{2});
                        elseif strcmp(tup{1}, 'inlink')
                            inlink = [];
                            for link = 2:length(tup)
                                inlink = [inlink; str2double(tup{link})];
                            end
                        elseif strcmp(tup{1}, 'outlink')
                            outlink = [];
                            for link = 2:length(tup)
                                outlink = [outlink; str2double(tup{link})];
                            end  
                        elseif strcmp(tup{1}, 'assignment_ratio')
                            ratio = [];
                            for link = 2:length(tup)
                                ratio = [ratio; str2double(tup{link})];
                            end  
                        elseif strcmp(tup{1}, 'junc_type')
                            junc_type = tup{2};
                        end
                    end
                    
                    % the last argument is T_grid. Not set yet, this field
                    % will be updated in the optimization program.
                    self.net.addJunc(junc_id, inlink, outlink, junc_type,...
                                     ratio, []);
                end
                
            end
            
            % delete(self.com.net_write_done);
            disp('Done read AIMSUN network information!');
      
        end
        
        
        
        %===============================================================
        % read historical boundary flow data
        % input:
        %       % data_file, the historical data file name
        %       % each row: 
        %           link_id,up/downstream,time (s),flow (veh/hr),speed
        %           (miles/hr)
        % output: saved in property historical_data
        %       % struct; his_data.link_1 fields: t, q_us, q_ds,
        %         v_us, v_ds
        function readHistoricalData(self, data_file)
            
            his_data = struct;
            
            content = self.readTextFile(data_file);
            
            % parse each line into numerical values
            for line = content
                
                % items = [link_id, link_bound, time, flow, speed]
                items = strsplit(line{1}, ',');
                
                linkStr = sprintf('link_%d', str2double(items{1}));
                
                if ~any(self.net.link_labels == str2double(items{1}))
                    error('ERROR: link %d not found in the network.\n',...
                        str2double(items{1}))
                end
                
                
                % initialize link field
                if ~isfield(his_data,linkStr)
                    his_data.(linkStr) = struct;
                end
                
                if strcmp(items{2}, 'upstream')
                    
                    if isfield(his_data.(linkStr), 't_us')
                        len = length(his_data.(linkStr).t_us);
                    else
                        % field has not been created yet
                        his_data.(linkStr).t_us=[];
                        his_data.(linkStr).q_us=[];
                        his_data.(linkStr).v_us=[];
                        len = 0;
                    end
                    his_data.(linkStr).t_us(len+1,1) = str2double(items{3});
                    his_data.(linkStr).q_us(len+1,1) = str2double(items{4})/3600;
                    
                    % the speed is optional
                    if length(items) == 5
                        his_data.(linkStr).v_us(len+1,1) = str2double(items{5})*1609/3600;
                    else
                        his_data.(linkStr).v_us(len+1,1) = NaN;
                    end
                    
                elseif strcmp(items{2}, 'downstream')
                    
                    if isfield(his_data.(linkStr), 't_ds')
                        len = length(his_data.(linkStr).t_ds);
                    else
                        % field has not been created yet
                        his_data.(linkStr).t_ds=[];
                        his_data.(linkStr).q_ds=[];
                        his_data.(linkStr).v_ds=[];
                        len = 0;
                    end
                    his_data.(linkStr).t_ds(len+1,1) = str2double(items{3});
                    his_data.(linkStr).q_ds(len+1,1) = str2double(items{4})/3600;
                    
                    % the speed is optional
                    if length(items) == 5
                        his_data.(linkStr).v_ds(len+1,1) = str2double(items{5})*1609/3600;
                    else
                        his_data.(linkStr).v_ds(len+1,1) = NaN;
                    end
                                        
                end
                
            end
            
            % save into property
            self.historical_data = his_data;
            
        end
        
        
        
        %===============================================================
        % warm up period
        % This period is for warming up. The reason is that the initial
        % density is normally not available or very unaccurate. Wait for
        % sometime until we got enough boundary flow data, then do MPC.
        % This function basically call getNewData to keep reading new data
        % generated by AIMSUN. It steps out once the run time at AIMSUN is
        % greater than warming_up_time.
        % input: 
        %       warming_up_time: the desired length of warming up time
        function warmUp(self, warming_up_time)
            
            self.dt_warm_up = warming_up_time;
            
            % wait until the warm up period has passed
            while self.t_now < self.dt_warm_up
                
                % Wating for the new data
                while ~self.getNewData();
                    pause(0.1); % pause 0.1 second
                end
                
                % got the new data; all data read and set in getNewData()
                % so do nothing here.
                
            end
            
            
        end
        
        
        %===============================================================
        % Get new data. It 
        % 1. returns false if new data is not available;
        % 2. returns true if new data available
        % 3. if new data avaliable, it will update the database in the 
        %    properties measurement_data, past_period_data
        % 4. It will update self.t_now to the latest time stamp
        function TF = getNewData(self)
            
            if ~exist(self.com.data_write_done,'file')
                % no new data
                TF = false;
                
            else
                % otherwise we got new data;
                % parse them and save in corresponding properties
                delete(self.com.data_write_done);
                content = self.readTextFile(self.com.data);
                delete(self.com.data);
                
                % parse each line into numerical values
                for line = content
                    
                    % items = [link_id, link_bound, time, flow, speed]
                    items = strsplit(line{1}, ',');
                    
                    linkStr = sprintf('link_%d', str2double(items{1}));
                    
                    if ~any(self.net.link_labels == str2double(items{1}))
                        error('ERROR: link %d not found in the network.\n',...
                            str2double(items{1}))
                    end
                    
                    % save to all_measurement_data
                    % initialize the link field
                    if ~isfield(self.all_measurement_data, linkStr)
                        self.all_measurement_data.(linkStr) = struct;
                    end
                    
                    if strcmp(items{2}, 'upstream')
                        
                        if isfield(self.all_measurement_data.(linkStr), 't_us')
                            len = length(self.all_measurement_data.(linkStr).t_us);
                        else
                            % field has not been created yet
                            self.all_measurement_data.(linkStr).t_us=[];
                            self.all_measurement_data.(linkStr).q_us=[];
                            self.all_measurement_data.(linkStr).v_us=[];
                            len = 0;
                        end
                        
                        self.all_measurement_data.(linkStr).t_us(len+1,1) = str2double(items{3});                        
                        % The flow in data is veh/hr, convert to veh/s for
                        % internal use
                        self.all_measurement_data.(linkStr).q_us(len+1,1) = str2double(items{4})/3600;
                        
                        % speed is in km/hr, convert to m/s for internal
                        if length(items) == 5
                            self.all_measurement_data.(linkStr).v_us(len+1,1) = str2double(items{5})*1000/3600;
                        else
                            self.all_measurement_data.(linkStr).v_us(len+1,1) = NaN;
                        end
                            
                            
                    elseif strcmp(items{2}, 'downstream')
                        
                        if isfield(self.all_measurement_data.(linkStr), 't_ds')
                            len = length(self.all_measurement_data.(linkStr).t_ds);
                        else
                            % field has not been created yet
                            self.all_measurement_data.(linkStr).t_ds=[];
                            self.all_measurement_data.(linkStr).q_ds=[];
                            self.all_measurement_data.(linkStr).v_ds=[];
                            len = 0;
                        end
                        
                        self.all_measurement_data.(linkStr).t_ds(len+1,1) = str2double(items{3});                        
                        self.all_measurement_data.(linkStr).q_ds(len+1,1) = str2double(items{4})/3600;
                        
                        % speed is optional
                        if length(items) == 5
                            self.all_measurement_data.(linkStr).v_ds(len+1,1) = str2double(items{5})*1000/3600;
                        else
                            self.all_measurement_data.(linkStr).v_ds(len+1,1) = NaN;
                        end
                        
                    end
                    
                end
                
                
                
                % update the current time, the sim start and end time.
                self.t_now = max([self.all_measurement_data.(linkStr).t_us;...
                               self.all_measurement_data.(linkStr).t_ds]);
                           
                self.t_sim_start = max(0, self.t_now - self.dt_past);
                self.t_sim_end = min(self.t_now + self.dt_predict, self.t_horizon_end);
                           
                % update the past_period_data and predict_period_data
                for link = self.net.link_labels'
                    
                    linkStr = sprintf('link_%d', link);
                    
                    % save the upstream past period data
                    if isfield(self.all_measurement_data, linkStr) && ...
                       isfield(self.all_measurement_data.(linkStr), 't_us')
                        % if we have those measurement data, then save them
                        % in the past period data property
                        index_past = self.all_measurement_data.(linkStr).t_us > self.t_sim_start &...
                                     self.all_measurement_data.(linkStr).t_us <= self.t_now;
                        self.past_period_data.(linkStr).t_us =...
                            self.all_measurement_data.(linkStr).t_us(index_past);
                        self.past_period_data.(linkStr).q_us =...
                            self.all_measurement_data.(linkStr).q_us(index_past);
                        self.past_period_data.(linkStr).v_us =...
                            self.all_measurement_data.(linkStr).v_us(index_past);
                    end
                    
                    % save the downstream past period data
                    if isfield(self.all_measurement_data, linkStr) && ...
                       isfield(self.all_measurement_data.(linkStr), 't_ds')
                        
                        index_past = self.all_measurement_data.(linkStr).t_ds > self.t_sim_start &...
                                     self.all_measurement_data.(linkStr).t_ds <= self.t_now;
                        self.past_period_data.(linkStr).t_ds =...
                            self.all_measurement_data.(linkStr).t_ds(index_past);
                        self.past_period_data.(linkStr).q_ds =...
                            self.all_measurement_data.(linkStr).q_ds(index_past);
                        self.past_period_data.(linkStr).v_ds =...
                            self.all_measurement_data.(linkStr).v_ds(index_past);
                    end
                    
                    % save the upstream predict period data from the
                    % historical data set
                    if isfield(self.historical_data, linkStr) && ...
                       isfield(self.historical_data.(linkStr), 't_us')
                        % if we have those measurement data, then save them
                        % in the past period data property
                        index_predict = self.historical_data.(linkStr).t_us > self.t_now &...
                                        self.historical_data.(linkStr).t_us < self.t_sim_end;
                        last_index = find( self.historical_data.(linkStr).t_us >= self.t_sim_end, 1 );
                        % t_us(end) must be self.t_sim_end
                        self.predict_period_data.(linkStr).t_us =...
                            [self.historical_data.(linkStr).t_us(index_predict);...
                             self.t_sim_end];
                        self.predict_period_data.(linkStr).q_us =...
                            [self.historical_data.(linkStr).q_us(index_predict);...
                             self.historical_data.(linkStr).q_us(last_index)];
                        self.predict_period_data.(linkStr).v_us =...
                            [self.historical_data.(linkStr).v_us(index_predict);...
                             self.historical_data.(linkStr).v_us(last_index)];
                    end
                    
                    % save the downstream predict period data from the
                    % historical data set
                    if isfield(self.historical_data, linkStr) && ...
                       isfield(self.historical_data.(linkStr), 't_ds')
                        % if we have those measurement data, then save them
                        % in the past period data property
                        index_predict = self.historical_data.(linkStr).t_ds > self.t_now &...
                                        self.historical_data.(linkStr).t_ds <= self.t_sim_end;
                        last_index = find( self.historical_data.(linkStr).t_ds >= self.t_sim_end, 1 );
                        self.predict_period_data.(linkStr).t_ds =...
                            [self.historical_data.(linkStr).t_ds(index_predict);...
                             self.t_sim_end];
                        self.predict_period_data.(linkStr).q_ds =...
                            [self.historical_data.(linkStr).q_ds(index_predict);...
                             self.historical_data.(linkStr).q_ds(last_index)];
                        self.predict_period_data.(linkStr).v_ds =...
                            [self.historical_data.(linkStr).v_ds(index_predict);...
                             self.historical_data.(linkStr).v_ds(last_index)];
                    end
                    
                    
                    
                    
                end
                
                % set the marker
                TF = true;
                
            end
        end
        
        
        %===============================================================
        % This function extracts measurement data in period, and return in
        % a boundary data struct to set the boundary conditions.
        % input: 
        %   period_start: the absolute time in the entire time horizon
        %   period_end: the absolute time in the entire time horizon
        % output:
        %   boundary_data: struct, BC_us are NOT normalized
        %       .(linkStr).T_us, T_us_cum, BC_us, T_ds, T_ds_cum, BC_ds
        function boundary_data = extractBoundaryData(self, period_start, period_end)
            
            % Only returns the measurement data. If also need the
            % historical data to set the boundary condition, then use the
            % function updateBoundaryCondition instead
            if period_end > self.t_now
                Warning('Use updateBoundaryCondition to extract the boundary data from historical database\n');
                return
            end
            
            bound_str_list = {'us', 'ds'};
            
            for link  = self.net.link_labels'
                
                linkStr = sprintf('link_%d', link);
                % Two boundaries, us, or ds
                for i = 1:length(bound_str_list)
                    
                    boundStr = bound_str_list{i};
                    % strings needed
                    t_bound = strcat('t_', boundStr);
                    q_bound = strcat('q_', boundStr);   
                    BC_bound = strcat('BC_', boundStr);  % by tradition, name q as BC
                    v_bound = strcat('v_', boundStr);
                    T_bound = strcat('T_', boundStr);
                    T_bound_cum = strcat('T_', boundStr, '_cum');
                    
                    % if the data is availabel in the database
                    if isfield(self.all_measurement_data, linkStr) && ...
                       isfield(self.all_measurement_data.(linkStr), t_bound)
                        
                        index_data = self.all_measurement_data.(linkStr).(t_bound) > period_start &...
                                     self.all_measurement_data.(linkStr).(t_bound) <= period_end;
                        t_meas =[ period_start;...
                            self.all_measurement_data.(linkStr).(t_bound)(index_data)];
                        
                        % add those data to the boundary data
                        boundary_data.(linkStr).(T_bound) = (t_meas(2:end) - t_meas(1:end-1));
                        boundary_data.(linkStr).(T_bound_cum) =...
                            [0; cumsum(boundary_data.(linkStr).(T_bound))];
                        boundary_data.(linkStr).(BC_bound) = self.all_measurement_data.(linkStr).(q_bound)(index_data);
                        boundary_data.(linkStr).(v_bound) = self.all_measurement_data.(linkStr).(v_bound)(index_data);
                        
                    else
                        % evenly discretize the simulation period to 30s
                        T_sim = period_start:30:period_end;
                        T_sim = unique([T_sim'; period_end]);
                        
                        boundary_data.(linkStr).(T_bound) = T_sim(2:end)-T_sim(1:end-1);
                        boundary_data.(linkStr).(T_bound_sum) = ...
                            [0; cumsum(boundary_data.(linkStr).(T_bound) ) ];
                        
                        % set data to be empty
                        boundary_data.(linkStr).(BC_bound) = [];
                        boundary_data.(linkStr).(v_bound) = [];
                        
                    end
                        
                    
                end
                
            end
            
        end
        
        %===============================================================
        % once data is available, then for each link
        % 1. Measurement and historical both availabel, then update
        % 2. Measurement data available, but historical not available, then
        %    update measurment data, and evenly discretize predict period
        %    and set flow as NaN by 30 s
        % 3. Measurement data not available, but historical data is
        %    available, then set all as historical data
        % 4. Neither is available, evenly discretize past and predict
        %    period by 30 s, data set as []
        function updateBoundaryCondition(self)
            
            bound_str_list = {'us', 'ds'};
            
            for link = self.net.link_labels'
                
                linkStr = sprintf('link_%d', link);
                
                % Two boundaries, us, or ds
                for i = 1:length(bound_str_list)
                    
                    boundStr = bound_str_list{i};
                    % strings needed
                    t_bound = strcat('t_', boundStr);
                    q_bound = strcat('q_', boundStr);   
                    BC_bound = strcat('BC_', boundStr);  % by tradition, name q as BC
                    v_bound = strcat('v_', boundStr);
                    T_bound = strcat('T_', boundStr);
                    T_bound_cum = strcat('T_', boundStr, '_cum');
                    
                    
                    % If both available
                    if isfield(self.past_period_data, linkStr) && ...
                            isfield(self.past_period_data.(linkStr), t_bound) && ...
                            isfield(self.predict_period_data, linkStr) && ...
                            isfield(self.predict_period_data.(linkStr), t_bound)
                    
                        % meas_data during self.t_sim_start : self.t_now
                        % t_bound are the time stamps of measurement
                        t_meas = [self.t_sim_start;...
                            self.past_period_data.(linkStr).(t_bound)];
                        % compute the duration of each time step
                        T_meas = (t_meas(2:end) - t_meas(1:end-1));
                        
                        % Same for historical data
                        t_his = [self.t_now;...
                            self.predict_period_data.(linkStr).(t_bound)];
                        T_his = t_his(2:end) - t_his(1:end-1);
                        
                        self.net.network_hwy.(linkStr).(T_bound) = [T_meas;T_his];
                        
                        % the relative starting time of each roll is 0.
                        self.net.network_hwy.(linkStr).(T_bound_cum) = [0;...
                            cumsum(self.net.network_hwy.(linkStr).(T_bound))];
                        
                        % the boundary flow measurement
                        self.net.network_hwy.(linkStr).(BC_bound) = ...
                            [self.past_period_data.(linkStr).(q_bound);...
                             self.predict_period_data.(linkStr).(q_bound)];
                         
                        % the velocity measurement
                        self.net.network_hwy.(linkStr).(v_bound) = ...
                            [self.past_period_data.(linkStr).(v_bound);...
                             self.predict_period_data.(linkStr).(v_bound)];  
                    
                    % If historical data is not available; and measurement
                    % data is available
                    elseif isfield(self.past_period_data, linkStr) && ...
                            isfield(self.past_period_data.(linkStr), t_bound) && ...
                            (~isfield(self.predict_period_data, linkStr) || ...
                             ~isfield(self.predict_period_data.(linkStr), t_bound))
                        
                        % meas_data during self.t_sim_start : self.t_now
                        % t_bound are the time stamps of measurement
                        t_meas = [self.t_sim_start;...
                            self.past_period_data.(linkStr).(t_bound)];
                        T_meas = (t_meas(2:end) - t_meas(1:end-1));
                        
                        % evenly discretize the time for historical data
                        t_his = self.t_now: 30 : self.t_sim_end;
                        t_his = unique([t_his'; self.t_sim_end]);
                        T_his = t_his(2:end) - t_his(1:end-1);
                        
                        self.net.network_hwy.(linkStr).(T_bound) = [T_meas;T_his];
                        
                        % the relative starting time of each roll is 0.
                        self.net.network_hwy.(linkStr).(T_bound_cum) = [0;...
                            cumsum(self.net.network_hwy.(linkStr).(T_bound))];
                        
                        % the boundary flow measurement
                        self.net.network_hwy.(linkStr).(BC_bound) = ...
                            [self.past_period_data.(linkStr).(q_bound);...
                             ones(length(T_his),1)*NaN];
                         
                        % the velocity measurement
                        self.net.network_hwy.(linkStr).(v_bound) = ...
                            [self.past_period_data.(linkStr).(v_bound);...
                             ones(length(T_his),1)*NaN];  
                        
                    % If measurement data is not available and historical
                    % data is available
                    elseif  (~isfield(self.past_period_data, linkStr) || ...
                            ~isfield(self.past_period_data.(linkStr), t_bound)) && ...
                            (isfield(self.predict_period_data, linkStr) && ...
                             isfield(self.predict_period_data.(linkStr), t_bound))
                        
                        % get all data from historical data set
                        index_sim = self.historical_data.(linkStr).(t_bound) > self.t_sim_start &...
                                        self.historical_data.(linkStr).(t_bound) < self.t_sim_end;
                        
                        last_index = find(self.historical_data.(linkStr).(t_bound) >= self.t_sim_end,1);
                        
                        T_sim = [self.t_sim_start; ...
                                 self.historical_data.(linkStr).(t_bound)(index_sim);...
                                 self.t_sim_end];
                        
                        % the time discretization
                        self.net.network_hwy.(linkStr).(T_bound) = ...
                            T_sim(2:end) - T_sim(1:end-1);
                        self.net.network_hwy.(linkStr).(T_bound_cum) = [0;...
                            cumsum(self.net.netwoek_hwy.(linkStr).(T_bound))];
                        
                        % the boundary flow measurement  
                        self.net.network_hwy.(linkStr).(BC_bound) = ...
                            [self.historical_data.(linkStr).(q_bound)(index_sim);
                            self.historical_data.(linkStr).(q_bound)(last_index)];
                        
                        % the velocity measurement
                        self.net.network_hwy.(linkStr).(v_bound) = ...
                            [self.historical_data.(linkStr).(v_bound)(index_sim);
                            self.historical_data.(linkStr).(v_bound)(last_index)];
                       
                        
                    % If neither is available
                    elseif  (~isfield(self.past_period_data, linkStr) || ...
                            ~isfield(self.past_period_data.(linkStr), t_bound)) && ...
                            (~isfield(self.predict_period_data, linkStr) || ...
                             ~isfield(self.predict_period_data.(linkStr), t_bound))
                        
                        % Evenly discretize the simulation period by 30 s
                        T_sim = self.t_sim_start: 30 : self.t_sim_end;
                        T_sim = unique([T_sim'; self.t_sim_end]);
                         
                        % the time discretization
                        self.net.network_hwy.(linkStr).(T_bound) = ...
                            T_sim(2:end) - T_sim(1:end-1);
                        self.net.network_hwy.(linkStr).(T_bound_cum) = [0;...
                            cumsum(self.net.netwoek_hwy.(linkStr).(T_bound))];
                        
                        % the boundary flow measurement  
                        self.net.network_hwy.(linkStr).(BC_bound) = [];
                        
                        % the velocity measurement
                        self.net.network_hwy.(linkStr).(v_bound) = [];
                                                
                    end
   
                end
                
 
                
            end
            
            % also need to update the junction grid
            for junc = self.net.junc_labels'
                
                juncStr = sprintf('junc_%d', junc);
                inlink = self.net.network_junc.(juncStr).inlabel(1);
                linkStr = sprintf('link_%d', inlink);
                
                % copy the link discretization here
                self.net.network_junc.(juncStr).T = ...
                    self.net.network_hwy.(linkStr).T_ds;
                self.net.network_junc.(juncStr).T_cum = ...
                    self.net.network_hwy.(linkStr).T_ds_cum;
                
            end
            
            
        end
        
        
        
        %===============================================================
        % update the grid of junctions to decrease the discretization
        % error. T_grid need to be updated both in network_jun and
        % network_hwy
        % input:
        %       T_grid: struct; link_1, the duration of each step
        %           T_grid.link_1.T_us or T_grid.link_1.T_ds
        %           T_grid.junc_1.T
        function updateBoundaryGrid(self, T_grid)
            
            % if T_grid is [], simply return
            % The grid has been set by updateBoundaryConditions
            if isempty(T_grid)
                
                return
                
            end
                        
            % iterate each field and set the corresponding junc or link
            fields = fieldnames(T_grid);
            for i = 1:lentgh(fields)
                
                nameStr = fields{i};
                
                % if it is a link
                if strcmp(nameStr(1:4), 'link')
                    
                    % update upstream grid
                    if isfield(T_grid.(nameStr),'T_us')
                        self.net.network_hwy.(nameStr).T_us =...
                            T_grid.(nameStr).T_us; 
                        self.net.network_hwy.(nameStr).T_us_cum = ...
                            [0; cumsum(T_grid.(nameStr).T_us)];
                        
                        % Note the boundary data should be emtpy
                        if ~isempty(self.net.network_hwy.(nameStr).q_us)
                            error('ERROR: Can not regrid upstream bound of link %s', nameStr);
                        end
                        
                    end
                    
                    % update downstream grid
                    if isfield(T_grid.(nameStr),'T_ds')
                        self.net.network_hwy.(nameStr).T_ds =...
                            T_grid.(nameStr).T_ds; 
                        self.net.network_hwy.(nameStr).T_ds_cum = ...
                            [0; cumsum(T_grid.(nameStr).T_ds)];
                        
                        % Note the boundary data should be emtpy
                        if ~isempty(self.net.network_hwy.(nameStr).q_ds)
                            error('ERROR: Can not regrid downstream bound of link %s', nameStr);
                        end
                        
                    end
                
                % if it is a junc
                elseif strcmp(nameStr(1:4), 'junc')
                   
                    self.net.network_junc.(nameStr).T = T_grid.(nameStr).T;
                    self.net.network_junc.(nameStr).T_cum = ...
                        [0; cumsum(T_grid.(nameStr).T)];
                end
                
            end
            
        end
        
        
        
        %===============================================================
        % update the initial condition of the links, in net property.
        % 1. At t = 0, IC can be left [], so no data constraints, this
        %    code will automatically discretize the link evenly to cells
        %    around 200 m
        % 2. Otherwise, IC need to be extracted from density estimation,
        %    which is implemented in postSolution class. NOTE, the IC is
        %    normalized value to kc on each link.
        % input:
        %       init_condition: struct, .(linkStr).X_grid_cum, float column
        %                               .(linkStr).IC, float column,
        function updateInitialCondition(self, init_condition)
            
            % if is nan, then evenly discretize each link into cells with
            % length around 200 m. Data will be set as NaN.
            if isempty(init_condition)
                
                for link = self.net.link_labels'
                    
                    linkStr = sprintf('link_%d', link);
                    
                    % the number of segments of around 200 m length
                    num_seg = ceil(self.net.network_hwy.(linkStr).para_postkm*1000/200);
                    
                    tmp_dx = self.net.network_hwy.(linkStr).para_postkm*1000/num_seg;
                    self.net.network_hwy.(linkStr).X_grid_cum = (0:num_seg)'*tmp_dx;
                    
                    % set the density data as NaN
                    self.net.network_hwy.(linkStr).IC = ones(num_seg,1)*NaN;
                    
                end
                
                return
                
            end
            
            % Otherwise, update the initial conditions
            fields = fieldnames(init_condition);
            
            % warning: the initial condition of some links is not provided
            if length(fields) ~= self.net.num_links
                warning('WARNING: the initial condition is not completely updated.\n')
            end
            
            % update each link
            for i = 1:length(fields)
                
                linkStr = fields{i};
                
                self.net.network_hwy.(linkStr).IC =...
                    init_condition.(linkStr).IC;
                
                % if grid provided, then set; otherwise not set,
                % initNetwork will discretize the space evenly.
                if isfield(init_condition.(linkStr), 'X_grid_cum')
                    self.net.network_hwy.(linkStr).X_grid_cum = ...
                        init_condition.(linkStr).X_grid_cum;
                end
                
            end
            
            
        end
       
        
        
        
        %===============================================================
        % This function applies the optimal control which is the optimal
        % continuous boundary flow
        % 1. AIMSUN meter takes flow parameters as the control input
        %    veh/hr. It actually approximate the flow using green red light.
        % 2. We need to make sure the flow is approximated right before
        %    changing the flow. E.g. 72 veh/hr is approximated by 49 s red
        %    and 1 s flash green. 
        % 3. It takes account of the realtime computation time. Only the 
        %    the control after t_now + dt_computation is applied  
        % 4. It saves the control signal in the format for AIMSUN to read
        %    start_time (s), duration (s), flow (veh/hr)
        % 5. AIMSUN is preset the flow in [0,900] veh/hr
        % input:
        %       x: the CP solution
        %       dt_computation: the computation time before applying the
        %           control
        %       dv_index: the decision variable index for exracting flow
        % output: the all_signal and signal_to_apply property will be
        %       updated
        %         the signal_to_apply will be written in file
        function applyControlByFlow(self, x, dt_computation, dv_index)
                        
            % for each link, update control signal
            for link = self.net.link_labels'
                
                linkStr = sprintf('link_%d', link);
                
                % if it is onramp, then control downstream flow
                if strcmp(self.net.network_hwy.(linkStr).para_linktype, 'onramp')
                
                    T_cum = self.net.network_hwy.(linkStr).T_ds_cum +...
                            self.t_sim_start;
                    
                    % round up the time.
                    t_signal_start = self.t_now + dt_computation;
                    
                    % find the continous flow from x                
                    index_flow = T_cum(2:end) > t_signal_start;
                    tmp_flow = x(dv_index.(linkStr).downstream(1,1):...
                                 dv_index.(linkStr).downstream(2,1));
                    % in veh/s
                    signal_flow = tmp_flow(index_flow);
                    
                    % find the setting time of each signal flow
                    tmp_setting_time = [t_signal_start;
                                           T_cum(T_cum > t_signal_start)];
                    signal_setting_time = tmp_setting_time(1:end-1);
                    signal_duration = tmp_setting_time(2:end) - tmp_setting_time(1:end-1);
                                
                    % discretize the signal into cycles
                    self.signal_to_apply = [signal_setting_time, ...
                                            signal_duration,...
                                            signal_flow];                
                    
                    % update all_signal property
                    if isempty(self.all_signal)
                        self.all_signal = self.signal_to_apply;
                    else
                        self.all_signal( self.all_signal(:,1) >= t_signal_start, :) = [];
                        self.all_signal = [self.all_signal; self.signal_to_apply];
                    end
                                            
                elseif strcmp(self.net.network_hwy.(linkStr).para_linktype, 'offramp')  
                    
                    % we have not defined an offramp actuator
                    error('ERROR: Offramp actuator not defined yet.\n')
                    
                end
                                

            end
            

            if ~isempty(self.signal_to_apply)
                fprintf('Time %d: Writing signals ...\n', self.t_now)
                
                fileID = fopen(self.com.signal,'w');
                
                % write header
                fprintf(fileID,'#signal_setting_time (s),cycle_duration (s), flow (veh/hr)\n');
                
                % write signal
                % NOTE: fprintf writes each column of the matrix as a row in
                % file, hence transpose the matrix.
                % convert the flow from veh/s to veh/hr, then write
                
                signal_to_apply_vehhr = self.signal_to_apply;
                signal_to_apply_vehhr(:,3) = signal_to_apply_vehhr(:,3)*3600;
                fprintf(fileID,'%.2f,%.2f,%d\n',(signal_to_apply_vehhr)');
                
                fclose(fileID);
                
                % write the flag file
                fileID = fopen(self.com.signal_write_done,'w');
                fprintf(fileID,'write done');
                fclose(fileID);
                fprintf('Wrote signal: %d ~ %d\n', signal_to_apply_vehhr(1,1),...
                    signal_to_apply_vehhr(end,1)+signal_to_apply_vehhr(end,2))
                
                % Make sure AIMSUN read the signal, signal file will be
                % deleted once AIMSUN finishes reading it.
                tmp_counter = 0;
                while exist(self.com.signal,'file')
                    if tmp_counter == 10
                        disp('Wait for AIMSUN to read the signal...\n')
                        tmp_counter = 0;
                    end
                    tmp_counter = tmp_counter+1;
                    pause(0.1)    % wait in second
                end
                
            else
                self.stopControl();
                
            end
            
        end
        
        
        
        
        %===============================================================
        % This function applies the optimal control which is discretized to
        % cycles.
        % 1. AIMSUN meter takes flow parameters as the control input
        %    veh/hr. It actually approximate the flow using green red light.
        % 2. We need to make sure the flow is approximated right before
        %    changing the flow. E.g. 72 veh/hr is approximated by 49 s red
        %    and 1 s flash green. 
        % 3. It takes account of the realtime computation time. Only the 
        %    the control after t_now + dt_computation is applied  
        % 4. It saves the control signal in the format for AIMSUN to read
        %    start_time (s), duration (s), flow (veh/hr)
        % 5. AIMSUN is preset the flow in [0,900] veh/hr
        % input:
        %       x: the CP solution
        %       links: int column vector; the links we would like to
        %           control
        %       dt_cycle: the control signal cycle. a new signal is set at 
        %           cycle 
        %       dt_computation: the computation time before applying the
        %           control
        % output: the all_signal and signal_to_apply property will be
        %       updated
        %         the signal_to_apply will be written in file
        % DO NOT USE.
        function applyControlByCycle(self, links, x, dt_cycle, dt_computation)
                        
            % for each link, update control signal
            for link = links'
                
                linkStr = sprintf('link_%d', link);
                
                % if it is onramp, then control downstream flow
                if strcmp(self.net.network_hwy.(linkStr).para_linktype, 'onramp')
                
                    T_cum = self.net.network_hwy.(linkStr).T_ds_cum +...
                            self.t_now;
                    
                    % round up the time.
                    t_signal_start = ceil((self.t_now + dt_computation)/dt_cycle)*...
                                        dt_cycle;
                    
                    % find the continous flow from x                
                    index_flow = T_cum(2:end) > t_signal_start;
                    tmp_flow = x(self.net.dv_index.(linkStr).downstream(1,1):...
                                    self.net.dv_index.(linkStr).downstream(2,1));
                    signal_flow = tmp_flow(index_flow);
                    
                    % find the setting time of each signal flow
                    tmp_setting_time = [t_signal_start;
                                           T_cum(T_cum > t_signal_start)];
                    signal_setting_time = tmp_setting_time(1:end-1);
                                
                    % discretize the signal into cycles
                    self.signal_to_apply = self.discretizeSignal(signal_setting_time,...
                                                signal_flow, dt_cycle);                
                    
                    % update all_signal property
                    self.all_signal( self.all_signal(:,1) >= t_signal_start, :) = [];
                    self.all_signal = [self.all_signal; self.signal_to_apply];
                    
                                            
                elseif strcmp(self.net.network_hwy.(linkStr).para_linktype, 'offramp')  
                    
                    % we have not defined an offramp actuator
                    error('ERROR: Offramp actuator not defined yet.\n')
                    
                else
                    % if not onramp or offramp, then not controllable
                    warning('WARNING: Link %d is not onramp or offramp, hence could not apply control.\n', link)
                    continue
                end
                                
                
                
            end
            
            
            % write signal to file
            disp('Writing signals...\n')
            
            fileID = fopen(self.com_signal_file,'w');
            
            % write header
            fprintf(fileID,'#signal_setting_time,cycle_duration,signal_value(green-0,red-1)\n');
            
            % write signal
            % NOTE: fprintf writes each column of the matrix as a row in
            % file, hence transpose the matrix.
            fprintf(fileID,'%.2f,%.2f,%d\n',self.signal_to_aply');
            
            fclose(fileID);
            
            % write the flag file
            fileID = fopen(self.com_signal_write_done,'w');
            fprintf(fileID,'write done');
            fclose(fileID);
            disp('Finished writing signals.\n')
            
            
        end
        
        
        
        % Discretize the boundary flow into cycles. 
        % At each cycle, we specify a new flow 
        % input: 
        %       signal_setting_time: nx1 vector, the set time of signal
        %       signal_flow: nx1 vector the flow (veh/s)
        %       dt_cycle: the cycle length
        % DO NOT USE
        function discretizeSignal(self, signal_setting_time, signal_flow, ...
                                    dt_cycle)
            
           error('ERROR: discretizeSignal not implemented yet.\n')                     
                                
        end
        
        
        
        %===============================================================
        % This function is called after MATLAB initialization. 
        % It create a file flag to tell AIMSUN that MATLAB is ready to
        % receive and process data (not ready to control yet)s
        function initMATLAB(self)
            
            % tell AIMSUN that the initialization of MATLAB is done, start
            % simulation
            dlmwrite(self.com.matlab_init_done, [], ';');
            fID = fopen(self.com.matlab_init_done);
            fclose(fID);
            
        end
        
        
        %===============================================================
        % This function is called when MATLAB is ready to start control 
        function startControl(self)
            
            % tell AIMSUN that the initialization of MATLAB is done, start
            % simulation
            dlmwrite(self.com.start_control, [], ';');
            fID = fopen(self.com.start_control);
            fclose(fID);
            
        end
        
        
        %===============================================================
        % This function is called when end time is reached and stop control 
        function stopControl(self)
            
            % tell AIMSUN that the initialization of MATLAB is done, start
            % simulation
            dlmwrite(self.com.stop_control, [], ';');
            fID = fopen(self.com.stop_control);
            fclose(fID);
            
        end
        
        
        
        %===============================================================
        % This function checks if the simulation in AIMSUN is completed.
        % output:
        %       True/False
        function TF = simulationCompleted(self)
            
            if exist(self.com.simulation_completed,'file')
                TF = true;
                return
            else
                TF = false;
                return
            end
            
        end
        
        
        
        %===============================================================
        % Utility function
        % read text file line by line. Save each line as a string with new
        % line character removed. The entire file is saved as a cell of strings. 
        % The line started with comment character % will be ignored
        % input: 
        %       filename: the filename to be read
        % output: 
        %       conten: cell of strings with each line as a string
        function content = readTextFile(~, filename)
            
            content = {};
            lineCounter = 0;
            
            if (~exist(filename,'file'))
                error('File %s does not exist.\n', filename)
            end
            
            fid = fopen(filename,'r');
            
            tline = fgetl(fid);
            while ischar(tline)
                
                if strcmp(tline(1),'%') || strcmp(tline(1),'#')
                    % comment line, ignore and continue 
                    tline = fgetl(fid);
                    continue
                end
                
                % otherwise a line, save in cell
                lineCounter = lineCounter + 1;
                content{lineCounter} = tline;
                tline = fgetl(fid);
            end
            
            fclose(fid);
            
        end
        
        
        
        %===============================================================
        % Utility function
        % This function plots the signal control output and the real data
        % measurement
        % input: 
        %       link: the link id of the on ramp, where the meter is
        %           installed at its downstream 
        function compareSignalandData(self, link)
            
            linkStr = sprintf('link_%d', link);
            data = [self.all_measurement_data.(linkStr).t_ds,...
                    self.all_measurement_data.(linkStr).q_ds];
            % convert to veh/hr
            data(:,2) = data(:,2)*3600;
            % shift data since the timestamp is when the data is collected    
            data(:,1) = [0; data(1:end-1,1)];
            if data(end,1) ~= self.t_horizon_end
                last_data = [self.t_horizon_end, data(end,2)];
                data = [data; last_data];
            end
            
            signal = self.all_signal(:,[1,3]);
            signal(:,2) = signal(:,2)*3600;
            if signal(end,1) ~= self.t_horizon_end
                last_signal = [self.t_horizon_end, signal(end,2)];
                signal = [signal; last_signal];
            end
            
            figure
            stairs( data(:,1), data(:,2),'b', 'LineWidth', 2);
            hold on
            stairs(signal(:,1), signal(:,2), 'r--', 'LineWidth',2);
            plot([self.dt_warm_up, self.dt_warm_up], [0, 1200], 'k', 'LineWidth', 2);
            legend('data','signal','control starts')
            ylabel('flow (veh/hr)','FontSize', 16);
            xlabel('time (s)', 'FontSize', 16);
            title('Applied signal and its effect', 'FontSize', 20);
            xlim([self.t_horizon_start, self.t_horizon_end])
            
            
        end
        
        
        
        
        %===============================================================
        % Utility function
        % This function fills the past_period_data with all measurement
        % data, then call other functions to construct CP and Mos to plot
        % 1. It fills the upstream and downstream data for all links up to
        %    t_now
        % 2. It clears the predict_data to empty
        function replayHorizon(self)
            
            % This can be done simply by setting 
            % t_sim_start = t_horizon_start; t_sim_end = t_now
            self.t_sim_start = self.t_horizon_start;
            self.t_sim_end = self.t_now;
            
            for link = self.net.link_labels'
                
                linkStr = sprintf('link_%d', link);
                
                % save the upstream past period data
                if isfield(self.all_measurement_data, linkStr) && ...
                        isfield(self.all_measurement_data.(linkStr), 't_us')
                    % if we have those measurement data, then save them
                    % in the past period data property
                    index_past = self.all_measurement_data.(linkStr).t_us > self.t_sim_start &...
                        self.all_measurement_data.(linkStr).t_us <= self.t_now;
                    self.past_period_data.(linkStr).t_us =...
                        self.all_measurement_data.(linkStr).t_us(index_past);
                    self.past_period_data.(linkStr).q_us =...
                        self.all_measurement_data.(linkStr).q_us(index_past);
                    self.past_period_data.(linkStr).v_us =...
                        self.all_measurement_data.(linkStr).v_us(index_past);
                end
                
                % save the downstream past period data
                if isfield(self.all_measurement_data, linkStr) && ...
                        isfield(self.all_measurement_data.(linkStr), 't_ds')
                    
                    index_past = self.all_measurement_data.(linkStr).t_ds > self.t_sim_start &...
                        self.all_measurement_data.(linkStr).t_ds <= self.t_now;
                    self.past_period_data.(linkStr).t_ds =...
                        self.all_measurement_data.(linkStr).t_ds(index_past);
                    self.past_period_data.(linkStr).q_ds =...
                        self.all_measurement_data.(linkStr).q_ds(index_past);
                    self.past_period_data.(linkStr).v_ds =...
                        self.all_measurement_data.(linkStr).v_ds(index_past);
                end
                
                
                % Clear the predict period data
                self.predict_period_data = struct;
               
            end
            
            
        end
        
        
        
        %===============================================================
        % Utility function
        % write all signal file for reply in AIMSUN
        function writeAllSignalReplay(self)
            
            fileID = fopen(self.com.all_signal,'w');
            
            % write header
            fprintf(fileID,'#signal_setting_time (s),cycle_duration (s), flow (veh/hr)\n');
            
            % write signal
            % NOTE: fprintf writes each column of the matrix as a row in
            % file, hence transpose the matrix.
            % convert the flow from veh/s to veh/hr, then write
            signal_to_apply_vehhr = self.all_signal;
            signal_to_apply_vehhr(:,3) = signal_to_apply_vehhr(:,3)*3600;
            fprintf(fileID,'%.2f,%.2f,%d\n',(signal_to_apply_vehhr)');
            
            fclose(fileID);
            
            % write the flag file
            fileID = fopen(self.com.all_signal_write_done,'w');
            fprintf(fileID,'write done');
            fclose(fileID);
            disp('Finished writing signals.')
            
        end
        
        
        %===============================================================
        % Utility function
        % read all data and replay the forward simulation
        function readAllDataReplay(self)
            
            if ~exist(self.com.all_data, 'file')
                error('all_data.txt does not exist.')
            end
            
            content = self.readTextFile(self.com.all_data);
            
            % parse each line into numerical values
                for line = content
                    
                    % items = [link_id, link_bound, time, flow, speed]
                    items = strsplit(line{1}, ',');
                    
                    linkStr = sprintf('link_%d', str2double(items{1}));
                    
                    if ~any(self.net.link_labels == str2double(items{1}))
                        error('ERROR: link %d not found in the network.\n',...
                            str2double(items{1}))
                    end
                    
                    % save to all_measurement_data
                    % initialize the link field
                    if ~isfield(self.all_measurement_data, linkStr)
                        self.all_measurement_data.(linkStr) = struct;
                    end
                    
                    if strcmp(items{2}, 'upstream')
                        
                        if isfield(self.all_measurement_data.(linkStr), 't_us')
                            len = length(self.all_measurement_data.(linkStr).t_us);
                        else
                            % field has not been created yet
                            self.all_measurement_data.(linkStr).t_us=[];
                            self.all_measurement_data.(linkStr).q_us=[];
                            self.all_measurement_data.(linkStr).v_us=[];
                            len = 0;
                        end
                        
                        self.all_measurement_data.(linkStr).t_us(len+1,1) = str2double(items{3});                        
                        % The flow in data is veh/hr, convert to veh/s for
                        % internal use
                        self.all_measurement_data.(linkStr).q_us(len+1,1) = str2double(items{4})/3600;
                        
                        % speed is in km/hr, convert to m/s for internal
                        if length(items) == 5
                            self.all_measurement_data.(linkStr).v_us(len+1,1) = str2double(items{5})*1000/3600;
                        else
                            self.all_measurement_data.(linkStr).v_us(len+1,1) = NaN;
                        end
                            
                            
                    elseif strcmp(items{2}, 'downstream')
                        
                        if isfield(self.all_measurement_data.(linkStr), 't_ds')
                            len = length(self.all_measurement_data.(linkStr).t_ds);
                        else
                            % field has not been created yet
                            self.all_measurement_data.(linkStr).t_ds=[];
                            self.all_measurement_data.(linkStr).q_ds=[];
                            self.all_measurement_data.(linkStr).v_ds=[];
                            len = 0;
                        end
                        
                        self.all_measurement_data.(linkStr).t_ds(len+1,1) = str2double(items{3});                        
                        self.all_measurement_data.(linkStr).q_ds(len+1,1) = str2double(items{4})/3600;
                        
                        % speed is optional
                        if length(items) == 5
                            self.all_measurement_data.(linkStr).v_ds(len+1,1) = str2double(items{5})*1000/3600;
                        else
                            self.all_measurement_data.(linkStr).v_ds(len+1,1) = NaN;
                        end
                        
                    end
                    
                end
            
            % update the current time, the sim start and sim end.
            self.t_now = max([self.all_measurement_data.(linkStr).t_us;...
                               self.all_measurement_data.(linkStr).t_ds]);
            self.t_sim_start = self.t_horizon_start;
            self.t_sim_end = self.t_now;
            
            % update the past_period_data and predict_period_data
            for link = self.net.link_labels'
                
                linkStr = sprintf('link_%d', link);
                
                % save the upstream past period data
                if isfield(self.all_measurement_data, linkStr) && ...
                        isfield(self.all_measurement_data.(linkStr), 't_us')
                    % if we have those measurement data, then save them
                    % in the past period data property
                    index_past = self.all_measurement_data.(linkStr).t_us > self.t_sim_start &...
                        self.all_measurement_data.(linkStr).t_us <= self.t_now;
                    self.past_period_data.(linkStr).t_us =...
                        self.all_measurement_data.(linkStr).t_us(index_past);
                    self.past_period_data.(linkStr).q_us =...
                        self.all_measurement_data.(linkStr).q_us(index_past);
                    self.past_period_data.(linkStr).v_us =...
                        self.all_measurement_data.(linkStr).v_us(index_past);
                end
                
                % save the downstream past period data
                if isfield(self.all_measurement_data, linkStr) && ...
                        isfield(self.all_measurement_data.(linkStr), 't_ds')
                    
                    index_past = self.all_measurement_data.(linkStr).t_ds > self.t_sim_start &...
                        self.all_measurement_data.(linkStr).t_ds <= self.t_now;
                    self.past_period_data.(linkStr).t_ds =...
                        self.all_measurement_data.(linkStr).t_ds(index_past);
                    self.past_period_data.(linkStr).q_ds =...
                        self.all_measurement_data.(linkStr).q_ds(index_past);
                    self.past_period_data.(linkStr).v_ds =...
                        self.all_measurement_data.(linkStr).v_ds(index_past);
                end
                
                % save the upstream predict period data from the
                % historical data set
                if isfield(self.historical_data, linkStr) && ...
                        isfield(self.historical_data.(linkStr), 't_us')
                    % if we have those measurement data, then save them
                    % in the past period data property
                    index_predict = self.historical_data.(linkStr).t_us > self.t_now &...
                        self.historical_data.(linkStr).t_us < self.t_sim_end;
                    last_index = find( self.historical_data.(linkStr).t_us >= self.t_sim_end, 1 );
                    % t_us(end) must be self.t_sim_end
                    self.predict_period_data.(linkStr).t_us =...
                        [self.historical_data.(linkStr).t_us(index_predict);...
                        self.t_sim_end];
                    self.predict_period_data.(linkStr).q_us =...
                        [self.historical_data.(linkStr).q_us(index_predict);...
                        self.historical_data.(linkStr).q_us(last_index)];
                    self.predict_period_data.(linkStr).v_us =...
                        [self.historical_data.(linkStr).v_us(index_predict);...
                        self.historical_data.(linkStr).v_us(last_index)];
                end
                
                % save the downstream predict period data from the
                % historical data set
                if isfield(self.historical_data, linkStr) && ...
                        isfield(self.historical_data.(linkStr), 't_ds')
                    % if we have those measurement data, then save them
                    % in the past period data property
                    index_predict = self.historical_data.(linkStr).t_ds > self.t_now &...
                        self.historical_data.(linkStr).t_ds <= self.t_sim_end;
                    last_index = find( self.historical_data.(linkStr).t_ds >= self.t_sim_end, 1 );
                    self.predict_period_data.(linkStr).t_ds =...
                        [self.historical_data.(linkStr).t_ds(index_predict);...
                        self.t_sim_end];
                    self.predict_period_data.(linkStr).q_ds =...
                        [self.historical_data.(linkStr).q_ds(index_predict);...
                        self.historical_data.(linkStr).q_ds(last_index)];
                    self.predict_period_data.(linkStr).v_ds =...
                        [self.historical_data.(linkStr).v_ds(index_predict);...
                        self.historical_data.(linkStr).v_ds(last_index)];
                end
                
                
                
                
            end
            
            
        end
        
        
        
    end
    
    
end
































