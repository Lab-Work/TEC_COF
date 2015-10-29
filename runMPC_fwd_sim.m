% Yanning Li, 
% Oct 14, 2015

% This script simply reads in the measurement data and forward simulate the
% entire horizon. The goal is to check if the estimated states matches the
% real states in AIMSUN, as well as how much the model error and
% measurement error affects the result.

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

% read all measurment data
meter.readAllDataReplay();

% remove link_330.upstream data since q3 = q1+q2;
meter.past_period_data.link_330.q_us = NaN*ones(length(meter.past_period_data.link_330.t_us),1);

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
errors.e_his = 0.2; % historical data error
errors.e_est = 0.1; % estimated initial condition error
errors.e_meas_flow = 0.05;
workzone_capacity = 0.6;
max_meter_rate = 900;   % veh/hr

% onramp queue limit
hard_queue_limit = struct;
% hard_queue_limit.link_390 = 1000; % onramp
% hard_queue_limit.link_330 = 300;
% hard_queue_limit.link_329 = 100;
% or queue_limit.(linkStr) which sets the limit of queue on freeways

% soft queue limit
soft_queue_limit = struct;
soft_queue_limit.link_330 = 100;

% update the boundary and initial condition
meter.updateBoundaryCondition();
meter.updateInitialCondition(init_condition);

%=========================================
% build the CP; set the matrix, and solve
CP = optProgram;
% start time, end time, queue_limit .(onramp) = 2/3 length
CP.setConfig(meter.net, 0, meter.t_now-meter.t_horizon_start,...
    meter.t_now,...
    hard_queue_limit, soft_queue_limit);

% constraints
CP.setConstraints(errors);
% CP.setWorkzoneCapacity([330], workzone_capacity);
% CP.setOnrampMeterMaxRate([390], max_meter_rate);

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

title_str = sprintf('No control %d ~ %d', meter.t_horizon_start, meter.t_sim_end);
        Mos.plotJuncs('all', title_str);
        
toc



q_329_us = x(CP.dv_index.link_329.upstream(1):CP.dv_index.link_329.upstream(2))*3600;
q_329_ds = x(CP.dv_index.link_329.downstream(1):CP.dv_index.link_329.downstream(2))*3600;
q_330_us = x(CP.dv_index.link_330.upstream(1):CP.dv_index.link_330.upstream(2))*3600;
q_330_ds = x(CP.dv_index.link_330.downstream(1):CP.dv_index.link_330.downstream(2))*3600;
q_390_us = x(CP.dv_index.link_390.upstream(1):CP.dv_index.link_390.upstream(2))*3600;
q_390_ds = x(CP.dv_index.link_390.downstream(1):CP.dv_index.link_390.downstream(2))*3600;

d_329_us = meter.all_measurement_data.link_329.q_us*3600;
d_329_ds = meter.all_measurement_data.link_329.q_ds*3600;
d_330_us = meter.all_measurement_data.link_330.q_us*3600;
d_330_ds = meter.all_measurement_data.link_330.q_ds*3600;
d_390_us = meter.all_measurement_data.link_390.q_us*3600;
d_390_ds = meter.all_measurement_data.link_390.q_ds*3600;

% check merge mass conservation
% [d_329_ds+d_390_ds, d_330_us, d_329_ds+d_390_ds-d_330_us]
[q_329_ds+q_390_ds, q_330_us, q_329_ds+q_390_ds-q_330_us]

% compare solution and data
% [q_329_ds, d_329_ds, q_329_ds-d_329_ds]
% [q_390_ds, d_390_ds, q_390_ds-d_390_ds]
% [q_330_ds, d_330_ds, q_330_ds-d_330_ds]







