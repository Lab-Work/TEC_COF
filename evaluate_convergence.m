% This script analyzes the convergence and stability of the convex junction
% solver in a merge example.
% Yanning Li, Feb 15, 2016.

% In theory, as the intervals in the boundary grid go to zero, the
% piecewise constant boudnary flow solution should be the same as the true
% solutioin.


tic
clearvars -except dbg

profile on
%===============================================================
%=================== Macro setting and Measures ================
% A new grid at each iteration. Save all the info
iter = struct;
iter.num_intervals = []; % column, number of intervals in a fixed horizon
iter.abs_error = [];     % column, absolute error of solution for each iter
iter.l2_error = [];      % the L2 error of boundary flow solution
iter.tt_time = [];       % total time for solving each iteration, including the overhead on building the LP
iter.solve_time = [];    % time needed for solving LP
iter.max_weight = [];    % the maximum weight coefficient in f, check for stability.

% time horizon
t_horizon_start = 0;
t_horizon_end = 150;

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

% Define an example network: link 1 & 2 merge to link 3
% length of each link
len_link1 = 1;    %km
len_link2 = 1;
len_link3 = 1;

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
% a temporary grid
sim_steps = 5;
step_length = (t_horizon_end-t_horizon_start)/sim_steps;   % seconds
T_init_grid = ones(sim_steps,1)*step_length;
net.addJunc(1, [1, 2]', 3, 'merge', [1; 1], T_init_grid);

% set initial conditoins
% Initial traffic density is initialized constants with even discretization
% Normalized to rho_c
% rho_tmp = randn(sim_steps, 1) + 0.5;
% rho_tmp(rho_tmp <= 0.2) = 0.2;
% rho_tmp(rho_tmp >= 1.5) = 1.5;

% rho_tmp = randn(sim_steps, 1) + 0.5;
% rho_tmp(rho_tmp <= 0.2) = 0.2;
% rho_tmp(rho_tmp >= 1.5) = 1.5;

% rho_tmp = randn(sim_steps, 1) + 0.5;
% rho_tmp(rho_tmp <= 0.2) = 0.2;
% rho_tmp(rho_tmp >= 1.5) = 1.5;

Ini.link_1.IC = net.network_hwy.link_1.para_kc*[1, 1, 1, 0, 0]';
Ini.link_2.IC = net.network_hwy.link_2.para_kc*[1, 1, 1, 0, 0]';
Ini.link_3.IC = net.network_hwy.link_3.para_kc*[1, 1, 1, 1, 1]';

net.setInitialCon(Ini);

% set the boundary condition
% Boundary condition at the two entrances and the one exit
% Normalized to q_max
% generate a random upstream flow
% q_tmp = randn(sim_steps, 1) + 1;
% q_tmp(q_tmp <= 0.5) = 0.5;
% q_tmp(q_tmp >= 0.9) = 0.9;
% q1_us_data = q_tmp;
q1_us_data = net.network_hwy.link_1.para_qmax*[1, 0, 1, 0, 1]';
    
% q_tmp = randn(sim_steps, 1) + 1;
% q_tmp(q_tmp <= 0.5) = 0.5;
% q_tmp(q_tmp >= 0.9) = 0.9;
% q2_us_data = q_tmp;
q2_us_data = net.network_hwy.link_2.para_qmax*[1, 0, 1, 0, 1]';

% q_tmp = randn(sim_steps, 1) + 1;
% q_tmp(q_tmp <= 0.5) = 0.5;
% q_tmp(q_tmp >= 0.9) = 0.9;
% q3_ds_data = q_tmp;
q3_ds_data = net.network_hwy.link_3.para_qmax*[1, 0, 1, 0, 1]';
% q3_ds_data = [];

% Set all parameters intended for control as empty or disabled
hard_queue_limit = struct;
soft_queue_limit = struct;

%===============================================================
%================== Compute the true solution ==================
% configuration paramters
sim_steps = 5;
step_length = (t_horizon_end-t_horizon_start)/sim_steps;   % seconds

% number of vehicles that we would consider as entropic solution
entropyTolerance = 0.01;

% This is the while loop for computing the true solution.
getEntropy = false;
loopCounter = 0;
T_junc = T_init_grid;
while getEntropy == false && loopCounter <= 50
    
    loopCounter = loopCounter+1;
    
    %==============================
    % update the grid at the junction
    net.network_junc.junc_1.T = T_junc;
    net.network_junc.junc_1.T_cum = [0; cumsum(T_junc)];
    
    %==============================                        
    % update boundary conditions along with the new grid
    % setBoundaryCon(obj, link, q_in, q_out, T_in, T_out)
    net.setBoundaryConForLink(1, q1_us_data, [], T_init_grid, T_junc);
    net.setBoundaryConForLink(2, q2_us_data, [], T_init_grid, T_junc);
    net.setBoundaryConForLink(3, [], q3_ds_data, T_junc, T_init_grid);
    
    %==============================
    % define and solve optimization program
    LP = optProgram;
    LP.setConfig(net, t_horizon_start, t_horizon_end, t_horizon_end,...
                 hard_queue_limit, soft_queue_limit)
    
    LP.setConstraints(errors);
    
    %==============================
    % Add objective functions
    LP.addEntropy(1);
    % LP.maxUpflow([1 2 4]);
    LP.maxDownflow(3, 1);
    % LP.maxError(1);
    %==============================    
    toc
    
    tic
    %solve program
    [x, fval, exitflag, output] = LP.solveProgram;
    
    %==============================
    % Post process the data
    Mos = postSolution(x, net, LP.dv_index, LP.end_time, dx_res, dt_res,...
                      hard_queue_limit, soft_queue_limit);
    % Mos.estimateState();
    % Mos.plotJuncs('all', 'Forward simulation with NO shockwave points constraints')
    
    [getEntropy, steps] = Mos.checkEntropy(entropyTolerance);
    
    if getEntropy == false
        hold on
        T_junc_all = Mos.updateTimeDiscretization(steps);
        % updated T_junc
        T_junc = T_junc_all.junc_1.T;
    end

end

% plot true solution
Mos.estimateState();
Mos.plotJuncs('all','true solution');

% aggregate and save the true solution
% sol1_true.T_cum = [t1, t2,...]
sol1_true = struct;
sol2_true = struct;

q1_true = x( LP.dv_index.link_1.downstream(1,1):...
    LP.dv_index.link_1.downstream(2,1) );
q2_true = x( LP.dv_index.link_2.downstream(1,1):...
    LP.dv_index.link_2.downstream(2,1) );

T_cum_1 = net.network_hwy.link_1.T_ds_cum(2:end);
T_cum_2 = net.network_hwy.link_2.T_ds_cum(2:end);

tmp_group = Mos.groupSameElement(q1_true, inf);
sol1_true.T_cum = T_cum_1(tmp_group(:,2));
sol1_true.q = tmp_group(:,3);

tmp_group = Mos.groupSameElement(q2_true, inf);
sol2_true.T_cum = T_cum_2(tmp_group(:,1));
sol2_true.q = tmp_group(:,3);



%===============================================================
%================== Run solver using different grid ============
for num_steps = 1:150
    
    fprintf('\nIteration %d', num_steps)
    
    % start timer
    dt_0 = now;
    
    step_length = (t_horizon_end-t_horizon_start)/num_steps;

    %==============================
    % update the grid at the junction
    T_junc = ones(num_steps,1)*step_length;
    T_junc_cum = [0; cumsum(T_junc)];
    T_junc_cum(end) = t_horizon_end;    % otherwise numerical error
    
    net.network_junc.junc_1.T = T_junc;
    net.network_junc.junc_1.T_cum = T_junc_cum;
    
    %==============================                        
    % update boundary conditions along with the new grid
    % setBoundaryCon(obj, link, q_in, q_out, T_in, T_out)
    net.setBoundaryConForLink(1, q1_us_data, [], T_init_grid, T_junc);
    net.setBoundaryConForLink(2, q2_us_data, [], T_init_grid, T_junc);
    net.setBoundaryConForLink(3, [], q3_ds_data, T_junc, T_init_grid);
    
    %==============================
    % define and solve optimization program
    LP = optProgram;
    LP.setConfig(net, t_horizon_start, t_horizon_end, t_horizon_end,...
                 hard_queue_limit, soft_queue_limit)
    
    LP.setConstraints(errors);
    
    %==============================
    % Add objective functions
    LP.addEntropy(1);
    LP.maxDownflow(3, 1);
    
    % timer lap 1
    dt_1 = now;
    
    % solve program
    [x, fval, exitflag, output] = LP.solveProgram;
    
    dt_2 = now;
    
    %==============================
    % save the solution
    sol1 = struct;
    sol2 = struct;
    
    tmp_T_cum = T_junc_cum(2:end);
    
    sol1.T_cum = tmp_T_cum;
    sol1.q = x( LP.dv_index.link_1.downstream(1,1):...
                LP.dv_index.link_1.downstream(2,1) );
    
    sol2.T_cum = tmp_T_cum;
    sol2.q = x( LP.dv_index.link_2.downstream(1,1):...
                LP.dv_index.link_2.downstream(2,1) );
    
    % since the priority paramter is 1, it suffices to compute only one abs
    % error and multiply by 2
    iter.num_intervals = [iter.num_intervals; num_steps];
    iter.abs_error = [iter.abs_error; 2*Mos.getAbsError(sol1_true, sol1) ];
    iter.tt_time = [iter.tt_time; (dt_2 - dt_0)*86400];
    iter.solve_time = [iter.solve_time; (dt_2 - dt_1)*86400];
    iter.max_weight = [iter.max_weight; max( abs(LP.f) )];
    
    fprintf('\nIteration %d: tt_time:%f, solve_time: %f\n',...
        num_steps, iter.tt_time(end), iter.solve_time(end));
    
end


%===============================================================
%======= Plot the convergence and stability result  ============
figure
plot(iter.num_intervals, iter.abs_error, 'LineWidth', 2);
grid on
xlabel('number of grid points', 'FontSize', 20);
ylabel('Abs error', 'FontSize', 20);
title('Convergence', 'FontSize', 24);

figure
[ax, h1, h2] = plotyy(iter.num_intervals, iter.tt_time,...
                      iter.num_intervals, iter.solve_time);
set(h1, 'color','red', 'LineWidth',2);
set(h2, 'color', 'blue', 'LineWidth', 2);
legend([h1, h2],'Computation time', 'Solve time');
xlabel('number of grid points', 'FontSize', 20);
title('Computation time', 'FontSize', 24);

figure
plot(iter.num_intervals, iter.max_weight, 'LineWidth', 2);
grid on
xlabel('number of grid points', 'FontSize', 20);
ylabel('Max coefficient weight', 'FontSize', 20);
title('Stability', 'FontSize', 24);

toc


% profile viewer
% p = profile('info');
% profsave(p,'./result/matlab_profiler_results');





