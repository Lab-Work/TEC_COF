% Yanning Li, Sep 07, this script is for running the model predictive
% control for an onramp. 

% TODO: to be finished. 

tic
clearvars -except dbg


meter = rampController_MPC;

%===============================================================
% Parameters for the simulation
% configController(self, simulation_end_time,...
%                                   dt_past, dt_predict, dt_roll_step)
meter.configController(2*60*60, 60*60, 30*60, 5*60,'onrampControl');

% set up communication
meter.setUpCommunication();

initEnv;

%===============================================================
% read the network information
meter.readNetworkFile();

% read historical data
meter.readHistoricalData('historical_data.txt');

%===============================================================
% start the simulation
meter.startSimulation();
meter.warmUp(30*60);

%===============================================================
% now start a rolling time horizon simulation
% while new data coming in and t_now
while (true)
    
    if getNewData()
        %=========================================
        meter.updateBounaryCondition()
        % simply set initial condition NaN since we have no data.
        meter.updateInitialCondition([])
        
        getEntropy = false;
        loopCounter = 0;
        
        while getEntropy == false && loopCounter <= 50
            
            % update the boundary discretization grid
            meter.updateBoundaryCondition();
            
            % update the initial condition based on the estimation
            meter.updateInitialCondition();
            
            % update the junction grid in net.junc_#.T
            meter.updateBoundaryGrid();
            
            %=========================================
            % build the CP; set the matrix, and solve
            CP = optProgram_MPC;
            CP.setConfig('onrampControl', t_start, t_end, t_now);
            
            % constraints
            CP.setConstraints(net, T_cum, T_dis);
            
            % objective function
            CP.addEntropy(net, 2, T_dis);
            % addOnRampControl(obj, network, juncs, T_dis)
            CP.addOnRampEntropy(net, 2, T_dis);
            CP.maxNetExitFlow(net);
            
            % solve the CP
            [x, fval, exitflag, output] = CP.solveProgram;
            
            %=========================================
            % post processing, check entropy and update discretization
            % function obj=postSolutionMPC(x)
            % do not solver the states inside links
            meter.Mos = postSolutionMPC(x);
            
            % visualize the result, (computationally heavy)
            % Mos.plotJuncs('all');
            
            [getEntropy, steps] = Mos.checkEntropy(entropyTolerance);
            
            if getEntropy == false
                hold on
                T_junc_cum = Mos.updateTimeDiscretization(steps);
                
                % a quick hack, only works for one junction
                T_junc = T_junc_cum.junc_1(2:end)-T_junc_cum.junc_1(1:end-1);
            end
            
        end
        
        %=========================================
        % plot the entropic result if needed
        Mos.plotJuncs('all');
      
        %=========================================
        % extract and format signal
        % signal defined as 0-green; 1-red; -1 inactive
        meter.applyControl()
      
    end
    
end

toc



















