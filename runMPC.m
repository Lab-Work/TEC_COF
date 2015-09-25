% Yanning Li, Sep 07, this script is for running the model predictive
% control for an onramp. 

% TODO: to be finished. 

tic
clearvars -except dbg


meter = rampController;

%===============================================================
% Parameters for the simulation
% configController(self, simulation_end_time,...
%                                   dt_past, dt_predict, dt_roll_step)
t_horizon_start = 0;
t_horizon_end = 2*60*60;
dt_past = 30*60;
dt_predict = 30*60;
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
queue_limit.onramp = 2/3; % of total length
% or queue_limit.(linkStr) which sets the limit of queue on freeways

%===============================================================
% start the simulation
meter.startSimulation();
dt_warm_up = 25*60; % 25 min warm up. Then at 30 min, the sensor 
meter.warmUp(dt_warm_up);

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
        if meter.t_now - meter.t_sim_start <= dt_past
            % if no need to roll
            init_condition = [];
        else
            init_condition = Mos.extractDensity(self.t_now ...
                - self.t_previous_roll_now);            
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
            % We decided to have a hard constraints on the length of queues
            % on the onramp.
            % CP.maxNetExitFlow(net);
            
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
        % Mos.plotJuncs('all');
      
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
        t_previous_roll_now = self.t_now;   % update the t_now in previous roll
        
    end
    
end

toc



















