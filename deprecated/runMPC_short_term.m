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
t_horizon_end = 2*60*60;
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
init_condition = [];

% T_grid is the discretization at the junctions
% If set as []. The toolbox will automatically discretize it into 
% cells with 30 s length
T_junc_grid = [];

% The percent errors of the data
errors = struct;
errors.e_default = 0.2;
errors.e_his = 0.3; % historical data error
errors.e_est = 0.2; % estimated initial condition error
errors.e_meas_flow = 0.1;

% onramp queue limit
queue_limit = struct;
queue_limit.onramp = 1/3; % of total length
% or queue_limit.(linkStr) which sets the limit of queue on freeways

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
while (~exist(meter.com.stop_control, 'file'))
    
    if meter.getNewData()
       
        % set a timer for the computation time
        dt_computation_start = now;
        
        % extract the initial condition from previous roll solution which
        % has relative time starting from 0.
        % at time self.t_now - self.t_previous_roll_now
        if meter.t_sim_start == 0
            % if no need to roll
            init_condition = [];
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
            CP.setConfig(meter.net, 0, meter.t_sim_end-meter.t_sim_start,...
                         queue_limit);
            
            % constraints
            CP.setConstraints(errors);
            
            % objective function
            CP.addEntropy('all');
            
            % a meaningful objective here. 
            % maximize the onramp flow while not causing congestion in the
            % dowanstream bottleneck
            CP.maxOnrampFlow('all');
            
            % solve the CP
            [x, fval, exitflag, output] = CP.solveProgram;
            
            %=========================================
            % post processing, check entropy and update discretization
            % function obj=postSolutionMPC(x)
            % do not solver the states inside links
            Mos = postSolution(x, meter.net, CP.dv_index,...
                               meter.t_sim_end-meter.t_sim_start,...
                               meter.dx_res, meter.dt_res);
            Mos.estimateState();
            
            % visualize the result, (computationally heavy)
            % Mos.plotJuncs('all');
            
%             [getEntropy, steps] = Mos.checkEntropy(entropyTolerance);
%             
%             if getEntropy == false
%                 hold on
%                 T_junc_grid = Mos.updateTimeDiscretization(steps);
%                 
%             end
             
        end
        
        %=========================================
        % plot the entropic result if needed
        Mos.plotJuncs('all');
%         f = gcf;
%         fig_index = fig_index+1;
%         saveas(f, sprintf('%d', fig_index), 'png');
%         close(f);
      
        % Get the real time computation time
        dt_computation_end = now;
        % dt_computation = (dt_computation_end - dt_computation_start)*86400;
        
        % for debugging, assume a static 15 second computation time.
        dt_computation = 15;    
        
        fprintf('The computation time is %f', dt_computation);
        
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



















