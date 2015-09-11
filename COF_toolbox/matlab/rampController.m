% Yanning Li, sep 07, 2015
% This class is for set up a model predictive controller
% It can be used for onramp metering of offramp control. For now, we only
% consider the work zone control scenario with only one merge or one
% diverge.
% The control scheme is as follows::
% 1. Model predictive control scheme. E.g. use last 1 hour boundary data, 
%    predict the next 30 min, assuming we have the 30 min historical data.
% 2. Roll forward once new measurement is available. The rolling step can
%    be 5 min.
% 3. Only the first period of the 30 min optimal control output will be
%    applied on the ramp. Once new optimal control is availabel (updated by
%    new measurement), the new optimal will be applied to the ramp.
% 4. This class simulates a real-time control. The computation time of the
%    solver is taken into account. Only the optimal control output after
%    t_now + dt_computation time will be applied.



classdef rampController < handle
    
    properties
        
        % The network to be controlled
        net;
        
        % The optimization program to be solved
        CP;
        
        % The post optimization object
        Mos;
        
        % historical data
        historical_data;
        
        % past period data
        past_period_data;
        
        % predict period data (extracted from historical data)
        predict_period_data;
        
        % all measurement data
        all_measurement_data;
        
        % all discretized control signal.
        all_signal;
        
        % the short period signal to be write to file and applied in AIMSUN
        signal_to_apply;
        
        % the starting time of the entire simulation. This is the abusolute
        % zero time of the system.
        t_horizon_start;
        
        % the end time of the entire simulation. The simulation will stop
        % at this time.
        t_horizon_end;
        
        % rolling step duration. This parameter should be the interval
        % between two measurement. Set this parameter to prevent updating
        % the control signal too fast
        dt_roll_step;
        
        % the duration for warming up.
        dt_warm_up;
        
        % Use the past dt_past measurement data to predict the future.
        % This value must be larger than dt_warm_up
        dt_past;
        
        % compute the optimal control for dt_predict in the future
        dt_predict;
        
        % the current time in system
        t_now;
        
        % other paramters that can be leave as default
        % the entropic solution tolerance in veh
        entropy_tolerance;  
        
        % The time space diagram resolution for visualization
        dx_res;     % m
        dt_res;     % s
        
        % the type of the controller: 'none', 'onramp', 'offramp'
        control_type;
        
        % accuracy of different types of data
        % The accuracy of the historical data
        e_historical_flow; 
        
        % The accuracy of new boundary flow data
        e_flow_measurement;
        e_speed_measurment;
        
        % The estimation error of the initial density
        e_init_estimation;
        
        % The following are flag files for communication between MATLAB and
        % AIMSUN
        % the network file written by AIMSUN, containing the network
        % topology and properties
        com_net_write_done;
        com_net_file;
        
        % the new boundary data file written by AIMSUN
        com_data_write_done;
        com_data_file;
        
        % the updated signal file written by MATLAB
        com_signal_write_done;  % if exists, meaning signal updated
        com_signal_file;
        
        % the initialization of MATLAB is done
        com_matlab_init_done;
        
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
            
            self.e_historical_flow = 0.3;
            self.e_flow_data = 0.05;
            self.e_speed_data = 0.1;
            
            self.e_init_estimation = 0.1;
            
        end
        
        
        %===============================================================
        % configure the control parameters
        function configController(self, simulation_end_time,...
                                  dt_past, dt_predict, dt_roll_step,...
                                  control_type)
                           
            self.t_horizon_end = simulation_end_time;
            self.dt_past = dt_past;
            self.dt_predict = dt_predict;
            self.dt_roll_step = dt_roll_step;
            
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
                items = strsplit(line,',');
                
                if strcmp(items{1}, 'net')
                    self.com_net_file = items{2};
                elseif strcmp(items{1}, 'net_write_done')
                    self.com_net_write_done = items{2};
                elseif strcmp(items{1}, 'data')
                    self.com_data_file = items{2};
                elseif strcmp(items{1}, 'data_write_done')
                    self.com_data_write_done = items{2};
                elseif strcmp(items{1}, 'signal')
                    self.com_signal_file = items{2};
                elseif strcmp(items{1}, 'signal_write_done')
                    self.com_signal_write_done = items{2};
                elseif strcmp(items{1}, 'matlab_init_done')
                    self.com_matlab_init_done = items{2};
                end

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
            while (~exist(self.com_net_write_done,'file'))
                pause(0.1);     % wait 0.5 second
            end
            disp('Got AIMSUN network information...');
            
            % remove flag file
            % should be deleted, keep here for debugging purpose
            delete(self.com_net_write_done);
            
            disp('Done read AIMSUN network information!');
            
            % set up the network
            self.net = initNetwork_MPC();
            
            % parse and set up the network.
            content = self.readTextFile(self.com_net_file);
            for line = content
                
                items = strsplit(line, ';');
                
                % if it is a link
                if strcmp(items{1},'link')
                    
                    % build a parameter struct
                    tmp_para = struct;
                    
                    % parse each property
                    for i = 2:length(items)
                        
                        tup = strsplit(items{i});
                        
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
                        
                        % initialize the link in the net
                        self.net.addLink( link_id,tmp_para,...
                                          num_lanes,len_km,...
                                          link_type)
                        
                    end
                    
                elseif strcmp(items{1},'junc')    
                    % if it is a junction
                    % parse each property
                    for i = 2:length(items)
                        
                        tup = strsplit(items{i});
                        
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
      
        end
        
        
        
        
        %===============================================================
        % read historical boundary flow data
        % input:
        %       % data_file, the historical data file name
        % output: saved in property historical_data
        %       % struct; his_data.link_1 fields: t, q_us, q_ds,
        %         v_us, v_ds
        function readHistoricalData(self, data_file)
            
            his_data = struct;
            
            content = self.readTextFile(data_file);
            
            % parse each line into numerical values
            for line = content
                
                % items = [link_id, link_bound, time, flow, speed]
                items = strsplit(line, ',');
                
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
                    his_data.(linkStr).q_us(len+1,1) = str2double(items{4});
                    his_data.(linkStr).v_us(len+1,1) = str2double(items{5});
                    
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
                    his_data.(linkStr).q_ds(len+1,1) = str2double(items{4});
                    his_data.(linkStr).v_ds(len+1,1) = str2double(items{5});
                    
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
            
            self.t_warm_up = warming_up_time;
            
            % wait until the warm up period has passed
            while self.t_now <= self.t_warm_up
                
                % Wating for the new data
                while ~self.getNewData();
                    pause(0.1); % pause 0.1 second
                end
                
                % got the new data; all data read and set in getNewData()
                % so do nothign here.
                
            end
            
            
        end
        
        
        %===============================================================
        % Get new data. It 
        % 1. return false if new data is not available;
        % 2. return true if new data available
        % 3. if new data avaliable, it will update the database in the 
        %    properties measurement_data, past_period_data
        function TF = getNewData(self)
            
            if ~exist(self.com_data_write_done,'file')
                % no new data
                TF = false;
                
            else
                % otherwise we got new data;
                % parse them and save in corresponding properties
                delete(self.com_data_write_done);
                content = self.readTextFile(self.com_data_file);
                delete(self.com_data_file);
                
                % parse each line into numerical values
                for line = content
                    
                    % items = [link_id, link_bound, time, flow, speed]
                    items = strsplit(line, ',');
                    
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
                        self.all_measurement_data.(linkStr).q_us(len+1,1) = str2double(items{4});
                        self.all_measurement_data.(linkStr).v_us(len+1,1) = str2double(items{5});
                        
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
                        self.all_measurement_data.(linkStr).q_ds(len+1,1) = str2double(items{4});
                        self.all_measurement_data.(linkStr).v_ds(len+1,1) = str2double(items{5});
                        
                    end
                    
                end
                
                % the past period start end time
                self.t_now = max(self.all_measurement_data.(linkStr).t_us,...
                               self.all_measurement_data.(linkStr).t_ds);
                t_past_period_start = max(0, self.t_now - self.dt_past);
                t_predict_period_end = min(self.t_now + self.dt_predict, self.t_horizon_end);
                           
                % update the past_period_data and predict_period_data
                for link = self.net.link_labels'
                    
                    linkStr = sprintf('link_%d', link);
                    
                    % save the upstream past period data
                    if isfield(self.all_measurement_data, linkStr) && ...
                       isfield(self.all_measurement_data.(linkStr), 't_us')
                        % if we have those measurement data, then save them
                        % in the past period data property
                        index_past = self.all_measurement_data.(linkStr).t_us > t_past_period_start &...
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
                        
                        index_past = self.all_measurement_data.(linkStr).t_ds > t_past_period_start &...
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
                                        self.historical_data.(linkStr).t_us <= t_predict_period_end;
                        self.predict_period_data.(linkStr).t_us =...
                            self.historical_data.(linkStr).t_us(index_predict);
                        self.predict_period_data.(linkStr).q_us =...
                            self.historical_data.(linkStr).q_us(index_predict);
                        self.predict_period_data.(linkStr).v_us =...
                            self.historical_data.(linkStr).v_us(index_predict);
                    end
                    
                    % save the upstream predict period data from the
                    % historical data set
                    if isfield(self.historical_data, linkStr) && ...
                       isfield(self.historical_data.(linkStr), 't_ds')
                        % if we have those measurement data, then save them
                        % in the past period data property
                        index_predict = self.historical_data.(linkStr).t_ds > self.t_now &...
                                        self.historical_data.(linkStr).t_ds <= t_predict_period_end;
                        self.predict_period_data.(linkStr).t_ds =...
                            self.historical_data.(linkStr).t_ds(index_predict);
                        self.predict_period_data.(linkStr).q_ds =...
                            self.historical_data.(linkStr).q_ds(index_predict);
                        self.predict_period_data.(linkStr).v_ds =...
                            self.historical_data.(linkStr).v_ds(index_predict);
                    end
                    
                    
                    
                    
                end
                
                % set the marker
                TF = true;
                
            end
        end
        
        
        %===============================================================
        % once data is available, then 
        % 1. set the boundary grid as t_now-dt_past: t_now+ dt_predict
        % 2. set the boundary condition in past_period by measurement data
        % 3. set the boundary conditino in predict_period by historical
        %    data
        function updateBoundaryCondition(self)
            
            t_past_period_start = max(0, self.t_now - self.dt_past);
            
            for link = self.net.link_labels'
                
                linkStr = sprintf('link_%d', link);
                % set the upstream
                if isfield(self.past_period_data, linkStr) && ...
                   isfield(self.past_period_data.(linkStr), 't_us')
                    
                    % t_us are the time stamps of measurement
                    t_data = [t_past_period_start;...
                              self.past_period_data.(linkStr).t_us];
                    tmp_T_us = (t_data(2:end) - t_data(1:end-1));
                   
                    if ~isfield(self.historical_data, linkStr) || ...
                       ~isfield(self.historical_data.(linkStr), 't_us')     
                        error('ERROR: the historical data for link %d is not set.\n',...
                            link)
                    end
                    
                    % compute the duration of each time step
                    t_data = [self.t_now;...
                              self.historical_data.(linkStr).t_us];
                    
                    self.net.network_hwy.(linkStr).T_us = [tmp_T_us;...
                        (t_data(2:end) - t_data(1:end-1))];
                    
                    % the relative starting time of each roll is 0.
                    self.net.network_hwy.(linkStr).T_us_cum = [0;...
                        cumsum(self.net.network_hwy.(linkStr).T_us)];
                    
                    % the boundary flow measurement
                    tmp_q_us = self.past_period_data.(linkStr).q_us;
                   
                    if ~isfield(self.historical_data, linkStr) || ...
                       ~isfield(self.historical_data.(linkStr), 'q_us')     
                        error('ERROR: the historical data for link %d is not set.\n',...
                            link)
                    end
                    
                    self.net.network_hwy.(linkStr).q_us = [tmp_q_us;...
                         self.historical_data.(linkStr).q_us];
                     
                    % the velocity measurement
                    tmp_v_us = self.past_period_data.(linkStr).v_us;
                   
                    if ~isfield(self.historical_data, linkStr) || ...
                       ~isfield(self.historical_data.(linkStr), 'v_us')     
                        error('ERROR: the historical data for link %d is not set.\n',...
                            link)
                    end
                    
                    self.net.network_hwy.(linkStr).v_us = [tmp_v_us;...
                         self.historical_data.(linkStr).v_us];
                     
                else
                    % this boundary does not have detector
                    % simply set the data as []. The discretization grid
                    % will be T_junc which will be set by function
                    % updateBoundaryGrid()
                    self.net.network_hwy.(linkStr).T_us = [];
                    self.net.network_hwy.(linkStr).T_us_cum = [];
                    self.net.network_hwy.(linkStr).q_us = [];
                    self.net.network_hwy.(linkStr).v_us = [];
                end
                
                
                % set the downstream boundary condition
                if isfield(self.past_period_data, linkStr) && ...
                   isfield(self.past_period_data.(linkStr), 't_ds')
                    % boundary data is available
                    
                    % t_us are the time stamps of measurement
                    t_data = [t_past_period_start;...
                              self.past_period_data.(linkStr).t_ds];
                    tmp_T_ds = (t_data(2:end) - t_data(1:end-1));
                   
                    if ~isfield(self.historical_data, linkStr) || ...
                       ~isfield(self.historical_data.(linkStr), 't_ds')     
                        error('ERROR: the historical data for link %d is not set.\n',...
                            link)
                    end
                    
                    % compute the duration of each time step
                    t_data = [self.t_now;...
                              self.historical_data.(linkStr).t_ds];
                    
                    self.net.network_hwy.(linkStr).T_ds = [tmp_T_ds;...
                        (t_data(2:end) - t_data(1:end-1))];
                    
                    % the relative starting time of each roll is 0.
                    self.net.network_hwy.(linkStr).T_ds_cum = [0;...
                        cumsum(self.net.network_hwy.(linkStr).T_ds)];
                    
                    % the boundary flow measurement
                    tmp_q_ds = self.past_period_data.(linkStr).q_ds;
                   
                    if ~isfield(self.historical_data, linkStr) || ...
                       ~isfield(self.historical_data.(linkStr), 'q_ds')     
                        error('ERROR: the historical data for link %d is not set.\n',...
                            link)
                    end
                    
                    self.net.network_hwy.(linkStr).q_ds = [tmp_q_ds;...
                         self.historical_data.(linkStr).q_ds];
                     
                    % the velocity measurement
                    tmp_v_ds = self.past_period_data.(linkStr).v_ds;
                   
                    if ~isfield(self.historical_data, linkStr) || ...
                       ~isfield(self.historical_data.(linkStr), 'v_ds')     
                        error('ERROR: the historical data for link %d is not set.\n',...
                            link)
                    end
                    
                    self.net.network_hwy.(linkStr).v_ds = [tmp_v_ds;...
                         self.historical_data.(linkStr).v_ds];
                     
                else
                    % this boundary does not have detector
                    % simply set the data as []. The discretization grid
                    % will be T_junc which will be set by function
                    % updateBoundaryGrid()
                    self.net.network_hwy.(linkStr).T_ds = [];
                    self.net.network_hwy.(linkStr).T_ds_cum = [];
                    self.net.network_hwy.(linkStr).q_ds = [];
                    self.net.network_hwy.(linkStr).v_ds = [];
                end
                
                
            end
            
            
        end
        
        
        
        %===============================================================
        % update the grid of junctions to decrease the discretization
        % error. T_junc need to be updated both in network_jun and
        % network_hwy
        % input:
        %       T_grid: struct; link_1, the duration of each step
        %       T_grid.link_1.T_us or T_grid.link_1.T_ds
        %       T_grid.junc_1.T
        function updateBoundaryGrid(self, T_grid)
            
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
        % 1. At t = 0, IC can be left [], so no data constraints
        % 2. Otherwise, IC need to be extracted from density estimation,
        %    which is implemented in postSolution class. NOTE, the IC is
        %    normalized value to kc on each link.
        % input:
        %       init_condition: struct, .(linkStr).X_grid_cum, float column
        %                               .(linkStr).IC, float column,
        %                               normalized to kc
        function updateInitialCondition(self, init_condition)
            
            fields = fieldnames(init_condition);
            
            % warning: the initial condition of some links is not provided
            if length(fields) ~= self.net.num_links
                warning('WARNING: the initial condition is not completely updated.\n')
            end
            
            % update each link
            for i = 1:length(fields)
                
                linkStr = fields{i};
                
                self.net.network_hwy.(linkStr).IC =...
                    init_condition.(linkStr).IC*...
                    self.net.network_hwy.(linkStr).para_kc;
                
                % if grid provided, then set; otherwise not set,
                % initNetwork will discretize the space evenly.
                if isfield(init_condition.(linkStr), 'X_grid_cum')
                    self.net.network_hwy.(linkStr).X_grid_cum = ...
                        init_condition.(linkStr).X_grid_cum;
                end
                
            end
            
            
        end
       
        
        
        
        
        %===============================================================
        % This function applies the optimal control 
        % 1. It converts the continuous control output to discrete 
        % control signals (e.g. gree-red)
        % 2. It takes account of the realtime computation time. Only the 
        %    the control after t_now + dt_computation is applied  
        % 3. It saves the control signal in the format for AIMSUN to read
        % input:
        %       x: the CP solution
        %       links: int column vector; the links we would like to
        %           control
        %       dt_cycle: the desired duration of each control status. 
        %           Boundary flow will be discretized in to cycles.
        %       dt_computation: the computation time before applying the
        %           control
        % output: the all_signal and signal_to_apply property will be
        %       updated
        %         the signal_to_apply will be written in file
        function applyControl(self, links, x, dt_cycle, dt_computation)
                        
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
        
        
        
        % Discretize the boundary flow into cycles. This function need to
        % be implemented after we get an idea how the ramp meter works in
        % AIMSUN.
        function discretizeSignal(self, signal_setting_time, signal_flow, ...
                                    dt_cycle)
            
           error('ERROR: discretizeSignal not implemented yet.\n')                     
                                
        end
        
        
        
        %===============================================================
        % This function is called after MATLAB initialization. 
        % It create a file flag to tell AIMSUN that MATLAB is ready to run
        % the simulation.
        function startSimulation(self)
            
            % tell AIMSUN that the initialization of MATLAB is done, start
            % simulation
            dlmwrite(self.com_matlab_init_done, [], ';');
            fID = fopen(self.com_matlab_init_done);
            fclose(fID);
            
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
        function content = readTextFile(self, filename)
            
            content = {};
            lineCounter = 0;
            
            if (~exist(filename,'file'))
                error('File %s does not exist.\n', filename)
            end
            
            fid = fopen('test.txt');
            
            tline = fgetl(fid);
            while ischar(tline)
                
                if strcmp(tline(1),'%')
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
        
        
        
        
    end
    
    
end
































