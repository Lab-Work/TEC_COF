% This script illustrates how to use the toolbox for a forward simulation
% at a merge jucntion.

tic
clearvars -except dbg

profile on
%===============================================================
% configuration paramters
t_horizon_start = 0;
sim_steps = 5;
step_length = 30;   % seconds
t_horizon_end = sim_steps*step_length;

% number of vehicles that we would consider as entropic solution
entropyTolerance = 1;

% Percent error of the full range
% e.g. [q_meas-e_meas_flow*q_max q_meas+e_meas_flow*q_max]
%      [rho_est-e_est*kc, rho_est+e_est*kc]
errors = struct;
errors.e_default = 0.0;
errors.e_his = 0.0; % historical data error
errors.e_est = 0.0; % estimated initial condition error
errors.e_meas_flow = 0.0;

% set the resolution for solving the HJ PDE on each link
dx_res = 2; % meters
dt_res = 2; % seconds

%===============================================================
% Define an example network: link 1 & 2 merge to link 3
% length of each link
len_link1 = 1.1;    %km
len_link2 = 1.2;
len_link3 = 1.3;

% default_para road parameters
default_para = struct;
default_para.beta_off = 0.2;
default_para.vf = 65*1609/3600;  %65 miles/hr
default_para.w = -(7.5349)*1609/3600;     %m/s calibrated from corsim
%default_para.w = -(8.125)*1609/3600;     %m/s Here v/w = 8, just to resolve the discritization error
default_para.kc_pl = 34.6569/1609;          %veh/m calibrated from corsim
default_para.qmax_pl = default_para.kc_pl*default_para.vf;  %veh/s
default_para.qon_max = default_para.kc_pl*default_para.vf; 
default_para.qoff_max = default_para.kc_pl*default_para.vf;
default_para.km_pl = default_para.kc_pl*(default_para.w-default_para.vf)/default_para.w;
default_para.v_min = 0*1609/3600;
default_para.v_max = 65*1609/3600;

net = initNetwork;
net.addLink(1, default_para, 1, len_link1, 'freeway');
net.addLink(2, default_para, 1, len_link2, 'freeway');
% 1.5 lanes simulate a limited capacity while still using the default_para
net.addLink(3, default_para, 1.5, len_link3, 'freeway');   
 
% function addJunc(self, junc, inlabel, outlabel, type_junc, ratio, T)
T_init_grid = ones(sim_steps,1)*step_length;
net.addJunc(1, [1, 2]', 3, 'merge', [1; 1], T_init_grid);

%===============================================================
% set initial conditoins
% Initial traffic density is initialized constants with even discretization
% Normalized to rho_c
% rho_tmp = randn(sim_steps, 1) + 0.5;
% rho_tmp(rho_tmp <= 0.2) = 0.2;
% rho_tmp(rho_tmp >= 1.5) = 1.5;
% Ini_1.IC = rho_tmp;
Ini_1.IC = [1, 1, 1, 0, 0]';

% rho_tmp = randn(sim_steps, 1) + 0.5;
% rho_tmp(rho_tmp <= 0.2) = 0.2;
% rho_tmp(rho_tmp >= 1.5) = 1.5;
% Ini_2.IC = rho_tmp;
Ini_2.IC = [1, 1, 1, 0, 0]';

% rho_tmp = randn(sim_steps, 1) + 0.5;
% rho_tmp(rho_tmp <= 0.2) = 0.2;
% rho_tmp(rho_tmp >= 1.5) = 1.5;
% Ini_3.IC = rho_tmp;
Ini_3.IC = [1, 1, 1, 1, 1]';

net.setInitialCon(1,Ini_1);
net.setInitialCon(2,Ini_2);
net.setInitialCon(3,Ini_3);

%===============================================================
% set the boundary condition
% Boundary condition at the two entrances and the one exit
% Normalized to q_max
% generate a random upstream flow
% q_tmp = randn(sim_steps, 1) + 1;
% q_tmp(q_tmp <= 0.5) = 0.5;
% q_tmp(q_tmp >= 0.9) = 0.9;
% q1_us_data = q_tmp;
q1_us_data = [1, 0, 1, 1, 1]';
    
% q_tmp = randn(sim_steps, 1) + 1;
% q_tmp(q_tmp <= 0.5) = 0.5;
% q_tmp(q_tmp >= 0.9) = 0.9;
% q2_us_data = q_tmp;
q2_us_data = [1, 0, 1, 1, 1]';

% q_tmp = randn(sim_steps, 1) + 1;
% q_tmp(q_tmp <= 0.5) = 0.5;
% q_tmp(q_tmp >= 0.9) = 0.9;
% q3_ds_data = q_tmp;
q3_ds_data = [1,1,1,1,1]';

%===============================================================
% Set all parameters intended for control as empty or disabled
hard_queue_limit = struct;
soft_queue_limit = struct;

%===============================================================
% This is the while loop for iteratively regridding and eventually getting
% the entropy condition
% Flags for iteratively updating time discretization
getEntropy = false;
loopCounter = 0;
T_junc = T_init_grid;
while getEntropy == false && loopCounter <=50
    
    loopCounter = loopCounter+1;
    
    %===============================================================
    % update the grid at the junction
    net.network_junc.junc_1.T = T_junc;
    net.network_junc.junc_1.T_cum = [0; cumsum(T_junc)];
    
    %===============================================================                         
    % update boundary conditions along with the new grid
    % setBoundaryCon(obj, link, q_in, q_out, T_in, T_out)
    net.setBoundaryCon(1, q1_us_data, [], T_init_grid, T_junc);
    net.setBoundaryCon(2, q2_us_data, [], T_init_grid, T_junc);
    net.setBoundaryCon(3, [], q3_ds_data, T_junc, T_init_grid);
    
    %===============================================================
    % define and solve optimization program
    LP = optProgram;
    LP.setConfig(net, t_horizon_start, t_horizon_end, t_horizon_end,...
                 hard_queue_limit, soft_queue_limit)
    
    LP.setConstraints( errors);
    
    %===============================================================
    % Add objective functions
    LP.addEntropy(1);
    % LP.maxUpflow([1 2 4]);
    % LP.minDownflow(net, 1);
    % LP.maxDownflow(net, [1; 2]);
    % LP.maxError(1);
    %===============================================================    
    toc
    
    tic
    %solve program
    [x, fval, exitflag, output] = LP.solveProgram;
    
    %===============================================================
    % Post process the data
    Mos = postSolution(x, net, LP.dv_index, LP.end_time, dx_res, dt_res,...
                      hard_queue_limit, soft_queue_limit);
    Mos.estimateState();
    % Mos.plotLinks('all')
    Mos.plotJuncs('all', 'Forward simulation with NO shockwave points constraints')
    
    [getEntropy, steps] = Mos.checkEntropy(entropyTolerance);
    
    if getEntropy == false
        hold on
        T_junc_cum = Mos.updateTimeDiscretization(steps);
        % updated T_junc
        T_junc = T_junc_cum.junc_1.T;
    end

end



toc


% profile viewer
% p = profile('info');
% profsave(p,'./result/matlab_profiler_results');





