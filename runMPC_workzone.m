% Yanning Li, 
% Oct 14, 2015

% This is a short term MPC control example.
% 1. The past period is 10 min, with 1 min data aggregation interval
% 2. The predicted period is 10 min in the future, with grid 1 min
% 3. Data is passed back every 5 min (5 rows of data)s
% 4. Assume the errors are 0, just for testing if the entropy condition
%   works
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
dt_past = 10*60;
dt_predict = 10*60;
dt_roll_step = 5*60;

meter.configController(t_horizon_start, t_horizon_end, dt_past, dt_predict, ...
                       'onrampControl');

% set up communication
meter.setUpCommunication('E:\\AIMSUN_MATLAB_COM\\COM_CONFIG.txt');


%===============================================================
% read the network information
meter.readNetworkFile();

% read historical data
meter.readHistoricalData('E:\\AIMSUN_MATLAB_COM\\historical_data.txt');

% no data for the initial condition, hence set as empty
% The code will automatically evenly discretized the link
% into cells with length smaller than 200 m
% NOTE: if we know the initial condition, then
% Here init_condition.(linkStr).IC is normalized
% init_condition = [];
init_condition.link_390.IC = zeros(5,1);
init_condition.link_390.X_grid_cum = [0:...
    meter.net.network_hwy.link_390.para_postkm*1000/5: ...
    meter.net.network_hwy.link_390.para_postkm*1000]';

init_condition.link_330.IC = zeros(5,1);
init_condition.link_330.X_grid_cum = [0:...
    meter.net.network_hwy.link_330.para_postkm*1000/5: ...
    meter.net.network_hwy.link_330.para_postkm*1000]';

init_condition.link_329.IC = zeros(5,1);
init_condition.link_329.X_grid_cum = [0:...
    meter.net.network_hwy.link_329.para_postkm*1000/5:...
    meter.net.network_hwy.link_329.para_postkm*1000]';


% T_grid is the discretization at the junctions
% If set as []. The toolbox will automatically discretize it into 
% cells with 30 s length
T_junc_grid = [];

% Percent error of the full range
% e.g. [q_meas-e_meas_flow*q_max q_meas+e_meas_flow*q_max]
%      [rho_est-e_est*kc, rho_est+e_est*kc]
errors = struct;
errors.e_default = 0.2;
errors.e_his = 0.3; % historical data error
errors.e_est = 0.1; % estimated initial condition error
errors.e_meas_flow = 0.05;

% onramp queue limit
hard_queue_limit = struct;
hard_queue_limit.link_390 = 1000; % onramp
% hard_queue_limit.link_330 = 300;
% hard_queue_limit.link_329 = 100;
% or queue_limit.(linkStr) which sets the limit of queue on freeways

% soft queue limit
soft_queue_limit = struct;
soft_queue_limit.link_330 = 300;

%===============================================================
% start the simulation
meter.initMATLAB();
dt_warm_up = 5*60;
meter.warmUp(dt_warm_up);

meter.startControl();
fig_index = 0;
%===============================================================
% now start a rolling time horizon simulation
% while new data coming in and t_now
while (~exist(meter.com.stop_control, 'file') && ...
       ~exist(meter.com.simulation_completed, 'file'))
    
    if meter.getNewData()
       
        % set a timer for the computation time
        dt_computation_start = now;
        
        % extract the initial condition from previous roll solution which
        % has relative time starting from 0.
        % at time self.t_now - self.t_previous_roll_now
        if meter.t_sim_start == 0
            % if no need to roll
            % do not update
            % init_condition = [];
        else
            init_condition = Mos.extractDensity(meter.t_sim_start ...
                - t_previous_sim_start);            
        end
        
        % update the boundary discretization grid
        meter.updateBoundaryCondition();
         
        % update the initial condition based on the estimation
        meter.updateInitialCondition(init_condition);
        
        % Now construct a
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
            CP.setWorkzoneCapacity([330], 0.7);
            
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
                               meter.t_sim_end-meter.t_sim_start,...
                               meter.dx_res, meter.dt_res,...
                               hard_queue_limit, soft_queue_limit);
            Mos.estimateState();
            
            % visualize the result, (computationally heavy)
            % Mos.plotJuncs('all');
            
%             [getEntropy, steps] = Mos.checkEntropy(entropyTolerance);
%             
%             if getEntropy == false
%                 hold on
%                 T_junc_grid = Mos.updateTimeDiscretization(steps);
%             end
             
        end
        
        %=========================================
        % plot the entropic result if needed
        title_str = sprintf('%d ~ %d', meter.t_sim_start, meter.t_sim_end);
        Mos.plotJuncs('all', title_str);
%         f = gcf;
%         fig_index = fig_index+1;
%         saveas(f, sprintf('%d', fig_index), 'png');
%         close(f);
      
        % Get the real time computation time
        dt_computation_end = now;
        % dt_computation = (dt_computation_end - dt_computation_start)*86400;
        
        % for debugging, assume a static 15 second computation time.
        dt_computation = 15;    
        
        fprintf('The computation time is %f\n', dt_computation);
        
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

toc



















