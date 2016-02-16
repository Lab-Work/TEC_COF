% Yanning Li, 
% Jan 10, 2016

% This is work zone onramp metering control example using the standard MPC
% link 329 (fwy) and 390 (onramp) merge to link 330 (fwy)
% 1. No past period is used. 
% 2. The predicted period is 10 min in the future, with grid 30 sec
% 3. Data is streaming back every 1 min.
% 4. Predicted horizon rolls forward 1 min each time when new data is
%   available.
% 5. Use the measurement data to estimate the current traffic states. 
% 6. Assume the errors are 0. Refer to our old paper for robust control

% Control objective:
% - Queue limit right before the downstream work zone
% - maximize on ramp flow with lower weight compared to the main line.

tic
clearvars -except dbg

meter = rampController;

%===============================================================
% Parameters for the simulation
% configController(self, simulation_end_time,...
%                                   dt_past, dt_predict, dt_roll_step)
t_horizon_start = 0;
t_horizon_end = 1*60*60;
dt_past = 0;    % no past period
dt_predict = 10*60;

meter.configController(t_horizon_start, t_horizon_end, dt_past, dt_predict, ...
                       'onrampControl');

% set up communication
meter.setUpCommunication('E:\\AIMSUN_MATLAB_COM\\COM_CONFIG.txt');


%===============================================================
% read the network information
meter.readNetworkFile();

% read historical data, assuming 
meter.readHistoricalData('E:\\AIMSUN_MATLAB_COM\\historical_data.txt');

% no data for the initial condition, hence set as empty
% The code will automatically evenly discretized the link
% into cells with length smaller than 200 m
% NOTE: if we know the initial condition, then
% Here init_condition.(linkStr).IC is normalized
% init_condition = [];
init_condition.link_390.IC = zeros(5,1);
init_condition.link_390.IC(5,1) = meter.net.network_hwy.link_390.para_kc*0.2;
init_condition.link_390.X_grid_cum = [0:...
    meter.net.network_hwy.link_390.para_postkm*1000/5: ...
    meter.net.network_hwy.link_390.para_postkm*1000]';

init_condition.link_330.IC = zeros(5,1);
init_condition.link_330.IC(5,1) = meter.net.network_hwy.link_330.para_kc*0.4;
init_condition.link_330.X_grid_cum = [0:...
    meter.net.network_hwy.link_330.para_postkm*1000/5: ...
    meter.net.network_hwy.link_330.para_postkm*1000]';

init_condition.link_329.IC = zeros(5,1);
init_condition.link_329.IC(5,1) = meter.net.network_hwy.link_329.para_kc*0.8;
init_condition.link_329.X_grid_cum = [0:...
    meter.net.network_hwy.link_329.para_postkm*1000/5:...
    meter.net.network_hwy.link_329.para_postkm*1000]';

% Make a copy of the net which includes the network profile and initial and
% boundary conditions for short-term estimation
net = meter.net;

% T_grid is the discretization at the junctions
% If set as []. The toolbox will automatically discretize it into 
% cells with 30 s length
T_junc_grid = [];

% Percent error of the full range
% e.g. [q_meas-e_meas_flow*q_max q_meas+e_meas_flow*q_max]
%      [rho_est-e_est*kc, rho_est+e_est*kc]
errors = struct;
errors.e_default = 0.0;
errors.e_his = 0.0; % historical data error
errors.e_est = 0.0; % estimated initial condition error
errors.e_meas_flow = 0.0;
workzone_capacity = 0.5;
max_meter_rate = 900;   % veh/hr

% hard queue limit
% exceeding this limit will cause infeasibility.
hard_queue_limit = struct;
% hard_queue_limit.link_390 = 1000; % onramp
% hard_queue_limit.link_330 = 300;
% hard_queue_limit.link_329 = 100;
% or queue_limit.(linkStr) which sets the limit of queue on freeways

% soft queue limit
% exceeding this limit will be penalized.
soft_queue_limit = struct;
soft_queue_limit.link_330 = 100;

%===============================================================
% start the simulation
meter.initMATLAB();
dt_warm_up = 0*60;
meter.warmUp(dt_warm_up);

t_previous_sim_start = 0;

meter.startControl();
%===============================================================
% now start a rolling time horizon simulation
% while new data coming in and t_now
while (~exist(meter.com.stop_control, 'file') && ...
       ~exist(meter.com.simulation_completed, 'file'))
    
    % if just finished the warm up or got new data
    % if meter.getNewData() || meter.t_now - t_horizon_start == dt_warm_up
    if meter.getNewData() 
        
        % entire horizon has been simulated. Stop control
        if meter.t_sim_start == t_horizon_end
            meter.stopControl();
            break;
        end
        
        % set a timer for the computation time
        dt_computation_start = now;
        
        % compute the initial condition
        if meter.t_sim_start ~= 0
            % solve a mini convex program to get the exact estimation of
            % the new initial condition for the predicted horizon
            LP = optProgram;
            tmp_soft_queue_limit = struct;
            tmp_hard_queue_limit = struct;
            
            % set the initial and boundary condition for the single step
            % estimation
            net.setInitialCon(init_condition);
            
            % extract the measurement data
            boundary_data = meter.extractBoundaryData(t_previous_sim_start,...
                meter.t_sim_start);
            % set the boundary condition of the mini CP
            net.setBoundaryCon(boundary_data);
            % Theoretically, q3 = q1+q2; however, measurement error exists.
            % Relax this to avoid infeasibility when error is set as 0
            net.network_hwy.link_330.BC_us = net.network_hwy.link_330.BC_us*NaN;
            % The junction boundary grid should also be upgraded
            net.network_junc.junc_428.T = net.network_hwy.link_330.T_us;
            net.network_junc.junc_428.T_cum = net.network_hwy.link_330.T_us_cum;
            
            LP.setConfig(net, 0, 0, ...
                meter.t_sim_start-t_previous_sim_start,...
                tmp_soft_queue_limit, tmp_hard_queue_limit);
            LP.setConstraints(errors);
            
            % This CP is only a layer over Berkeley HJ PDE Solver for easy 
            % implementation, and gives exact solution, maximization here
            % will not change the solution, just for consistency of CPLEX
            LP.maxDownflow([330], 1);
            % solve the CP
            [x, fval, exitflag, output] = LP.solveProgram;
            
            % post processing, check entropy and update discretization
            % function obj=postSolutionMPC(x)
            % do not solver the states inside links
            State = postSolution(x, net, LP.dv_index,...
                               meter.t_sim_start-t_previous_sim_start,...
                               2, 2,...
                               tmp_hard_queue_limit, tmp_soft_queue_limit);
%             State.estimateState();
%             State.plotJuncs('all', sprintf('State estimate for %d',...
%                                            meter.t_sim_start ));
            
            init_condition = State.extractDensity(meter.t_sim_start ...
                - t_previous_sim_start);            
        end
        
        % Now solve the predicted horizon control problem
        % Use historical data to update the boundary conditions 
        meter.updateBoundaryCondition();
         
        % update the initial condition based on the estimation
        meter.updateInitialCondition(init_condition);
        
        getEntropy = false;
        loopCounter = 0;
        
        % disable the grid updating part.
        while getEntropy == false && loopCounter <= 0
            
            loopCounter = loopCounter + 1;
            
            % update the junction grid in net.junc_#.T
            % This is for adaptive griding.
            meter.updateBoundaryGrid(T_junc_grid);
            
            %=========================================
            % build the CP; set the matrix, and solve
            CP = optProgram;
            % start time, end time, queue_limit .(onramp) = 2/3 length
            CP.setConfig(meter.net, 0, meter.t_now-meter.t_sim_start,...
                         meter.t_sim_end-meter.t_sim_start,...
                         hard_queue_limit, soft_queue_limit);
            
            % constraints
            CP.setConstraints(errors);
            CP.setWorkzoneCapacity([330], workzone_capacity);
            CP.setOnrampMeterMaxRate([390], max_meter_rate);
            
%             % f v1
%             % The order is important
%             CP.maxDownflow([330], 1);
%             CP.penalizeCongestion;
%             CP.applyEntropy;
%             CP.maxOnrampFlow('all');
            
            % f v2
            % The order is important
            CP.penalizeCongestion_v2(1);
            CP.maxOnrampFlow_v2(390, 1);
            CP.applyEntropy_v2;
            
            
            % solve the CP
            [x, fval, exitflag, output] = CP.solveProgram;
            
            %=========================================
            % post processing, check entropy and update discretization
            % function obj=postSolutionMPC(x)
            % do not solver the states inside links
%             Mos = postSolution(x, meter.net, CP.dv_index,...
%                                meter.t_sim_end-meter.t_sim_start,...
%                                meter.dx_res, meter.dt_res,...
%                                hard_queue_limit, soft_queue_limit);
%             Mos.estimateState();
%             
%             [getEntropy, steps] = Mos.checkEntropy(entropyTolerance);
%             
%             if getEntropy == false
%                 hold on
%                 T_junc_grid = Mos.updateTimeDiscretization(steps);
%             end
             
        end
        
        %=========================================
        % plot the entropic result if needed
%         title_str = sprintf('%d ~ %d', meter.t_sim_start, meter.t_sim_end);
%         Mos.plotJuncs('all', title_str);
      
        % Get the real time computation time
        dt_computation_end = now;
        realtime_dt_computation = (dt_computation_end - dt_computation_start)*86400;
        % dt_computation = (dt_computation_end - dt_computation_start)*86400;
        
        % for debugging, assume a static 15 second computation time.
        dt_computation = 15;    
        
        fprintf('The computation time is %f\n', realtime_dt_computation);
        
        %=========================================
        % extract and format signal
        % apply the control as a flow meter
        meter.applyControlByFlow(x, dt_computation, CP.dv_index);
        
        %=========================================
        % save the state solution which will be used for extracting the new
        % initial condition
        flow_sol = x;
        t_previous_roll_now = meter.t_now;   % update the t_now in previous roll
        t_previous_sim_start = meter.t_sim_start;
    end
    
end

% write all the signal into file, then reply in AIMSUN.
meter.writeAllSignalReplay;

%==========================================================================
% Plot the entire simulation horizon

% set the past period to be the entire time horizon instead of dt_past
meter.replayHorizon();

% update the boundary discretization grid
meter.past_period_data.link_330.q_us = NaN*ones(length(meter.past_period_data.link_330.t_us),1);

meter.updateBoundaryCondition();

% update the initial condition based on the estimation
init_condition.link_390.IC =  meter.net.network_hwy.link_390.para_kc*ones(5,1)*0.2;
init_condition.link_390.X_grid_cum = [0:...
    meter.net.network_hwy.link_390.para_postkm*1000/5: ...
    meter.net.network_hwy.link_390.para_postkm*1000]';

init_condition.link_330.IC =  meter.net.network_hwy.link_330.para_kc*ones(5,1)*0.2;
init_condition.link_330.X_grid_cum = [0:...
    meter.net.network_hwy.link_330.para_postkm*1000/5: ...
    meter.net.network_hwy.link_330.para_postkm*1000]';

init_condition.link_329.IC =  meter.net.network_hwy.link_329.para_kc*ones(5,1)*0.2;
init_condition.link_329.X_grid_cum = [0:...
    meter.net.network_hwy.link_329.para_postkm*1000/5:...
    meter.net.network_hwy.link_329.para_postkm*1000]';

meter.updateInitialCondition(init_condition);

%=========================================
% build the CP; set the matrix, and solve
CP = optProgram;
% start time, end time, queue_limit .(onramp) = 2/3 length
CP.setConfig(meter.net, 0, meter.t_now-meter.t_horizon_start,...
    meter.t_now,...
    hard_queue_limit, soft_queue_limit);

% constraints
% Here e_est is only used to constrain initial condition which we know is 0
% errors.e_est = 0.5;
% errors.e_meas_flow = 0.04;
CP.setConstraints(errors);
CP.setWorkzoneCapacity([330], workzone_capacity);
CP.setOnrampMeterMaxRate([390], max_meter_rate);

% a meaningful objective here.
% The order is important
CP.maxDownflow([330], 1);
CP.penalizeCongestion;
CP.applyEntropy;
CP.maxOnrampFlow('all');

% solve the CP
[x, fval, exitflag, output] = CP.solveProgram;

%=========================================
% post processing, check entropy and update discretization
% function obj=postSolutionMPC(x)
% do not solver the states inside links
Mos = postSolution(x, meter.net, CP.dv_index,...
    meter.t_sim_end-meter.t_horizon_start,...
    meter.dx_res, meter.dt_res,...
    hard_queue_limit, soft_queue_limit);
Mos.estimateState();

title_str = sprintf('%d ~ %d', meter.t_horizon_start, meter.t_sim_end);
        Mos.plotJuncs('all', title_str);

meter.compareSignalandData([390])
        
toc



















