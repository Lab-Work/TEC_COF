%% An example for solving a merge junction using a convex program
% This script illustrates how to use the TEC toolbox for solving an
% unsignalized merge junction.

%% Configure parameters

clearvars -except dbg

% profile on

t_horizon_start = 0;
sim_steps = 5;
step_length = 30;   % seconds
t_horizon_end = sim_steps*step_length;

%%
% A tolerance threshold for identifying the exact solution.
% If the difference between the obtained solution and the exact solution 
% in each time interval is less than 1 veh, then it is considered as the
% exact solution.
exactTolerance = 1;

%%
% Percent error of the full range. For example:
%
% $$[q_{meas}-e_{measflow}\times q_{max} q_{meas}+e_{measflow}\times
% q_{max}]$$
% 
% $$[\rho_{est}-e_{est}\times \rho_c, \rho_{est}+e_{est}\times \rho_c]$$
%
% In this merge solver example, the initial and boundary condition errors
% are assumed to be exact. Errors can be added to accommodate measurement noise
% in traffic estimation applications.
errors = struct;
errors.e_default = 0.0; % default error used if others are not defined.
errors.e_his = 0.0; % historical flow data error
errors.e_est = 0.0; % estimated initial condition error
errors.e_meas_flow = 0.0;   % measurement flow data error

%%
% Set the resolution for solving the HJ PDE on each link:
dx_res = 100; % meters
dt_res = 10; % seconds

%%
% Set the default_para road parameters:
default_para = struct;
default_para.beta_off = 0.2;
default_para.vf = 65*1609/3600;  %65 miles/hr
default_para.w = -(7.5349)*1609/3600;     %m/s calibrated from AIMSUN
default_para.kc_pl = 34.6569/1609;        %veh/m calibrated from AIMSUN
default_para.qmax_pl = default_para.kc_pl*default_para.vf;  %veh/s
default_para.qon_max = default_para.kc_pl*default_para.vf; 
default_para.qoff_max = default_para.kc_pl*default_para.vf;
default_para.km_pl = default_para.kc_pl*(default_para.w-default_para.vf)/default_para.w;
default_para.v_min = 0*1609/3600;
default_para.v_max = 80*1609/3600;

%% Define an example merge network: 
%%
% The link 1 & 2 merge to link 3, each with the link length in km:
len_link1 = 1.1;       % km
len_link2 = 1.2;
len_link3 = 1.3;

%% 
% Define the network as an initNetowrk object. Refer to the iniNetwork doc
% for details.
net = initNetwork;
net.addLink(1, default_para, 1, len_link1, 'freeway');
net.addLink(2, default_para, 1, len_link2, 'freeway');
% 1.5 lanes simulate a limited capacity while still using the default_para
net.addLink(3, default_para, 1.5, len_link3, 'freeway');   

T_init_grid = ones(sim_steps,1)*step_length;
net.addJunc(1, [1, 2]', 3, 'merge', [1; 1], T_init_grid);

%% Set the initial and boundary conditions
% Set initial conditoins randomly or staticly.
% Initial traffic density is initialized constants with even
% discretization.

Ini.link_1.IC = net.network_hwy.link_1.para_kc*[1, 1, 1, 0, 0]';
Ini.link_2.IC = net.network_hwy.link_2.para_kc*[1, 1, 1, 0, 0]';
Ini.link_3.IC = net.network_hwy.link_3.para_kc*[1, 1, 1, 1, 1]';

net.setInitialCon(Ini);

%%
% Set the boundary condition. The downstream boundary flow is not set to
% aviod infeasibility.
q1_us_data = net.network_hwy.link_1.para_qmax*[1, 0, 1, 1, 1]';
q2_us_data = net.network_hwy.link_2.para_qmax*[1, 0, 1, 1, 1]';
q3_ds_data = net.network_hwy.link_3.para_qmax*[1, 0, 1, 0, 1]';

%%
% Queue limits are used in traffic control applications. In this merge
% solver, simply set them as empty struct.
hard_queue_limit = struct;
soft_queue_limit = struct;

%% Solve the internal boundary flows exactly
% The proposed numerical scheme can compute the internal boundary flows
% using a singel convex program no a given time grid over the entire time 
% horizon. To obtain the exact solutions, an iterative binary search
% algorithm is used to locate the intersection of the shockwaves at the
% junction.
% Each loop computes the internal boundary flow solution; checks if it is
% the exact solution; and updates the time discretization grid at the
% junction.

gotExact = false;   % Flag for identifying if the exact solution is obtained
loopCounter = 0;    % A loop counter to stop the iteration.
T_junc = T_init_grid;   % Start with an initial tie grid.
while gotExact == false && loopCounter <=50
    
    loopCounter = loopCounter+1;
    
    %===============================================================
    % set the grid at the junction
    net.network_junc.junc_1.T = T_junc;
    net.network_junc.junc_1.T_cum = [0; cumsum(T_junc)];
    
    %===============================================================                         
    % update boundary conditions along with the new grid
    % setBoundaryCon(obj, link, q_in, q_out, T_in, T_out)
    net.setBoundaryConForLink(1, q1_us_data, [], T_init_grid, T_junc);
    net.setBoundaryConForLink(2, q2_us_data, [], T_init_grid, T_junc);
    net.setBoundaryConForLink(3, [], q3_ds_data, T_junc, T_init_grid);
    
    %===============================================================
    % define and solve optimization program
    LP = optProgram;
    LP.setConfig(net, t_horizon_start, t_horizon_end, t_horizon_end,...
                 hard_queue_limit, soft_queue_limit)
    
    LP.setConstraints( errors);
    
    %===============================================================
    % Add objective functions
    LP.applyAdmissibleCon(1);
    
    %===============================================================    
    % solve program
    [x, fval, exitflag, output] = LP.solveProgram;
    
    %===============================================================
    % post process the solution using postSolution class
    Mos = postSolution(x, net, LP.dv_index, LP.end_time, dx_res, dt_res,...
                      hard_queue_limit, soft_queue_limit);
                  
    % Uncommente the following two lines to plot the solution in each iteration.
    % Mos.estimateState();
    % Mos.plotJuncs('all', 'Solving a merge')
    
    % check if the solution is exact
    [gotExact, steps] = Mos.checkSolution(exactTolerance);
    
    % if the solution is not exact, update the time grid
    if gotExact == false
        T_junc_cum = Mos.updateTimeGrid(steps);
        T_junc = T_junc_cum.junc_1.T;
    end

end

%% Output and visualization
% Visualize the final exact solution. 
Mos.estimateState();
Mos.plotJuncs('all', 'Solve a merge')


% profile viewer
% p = profile('info');
% profsave(p,'./result/matlab_profiler_results');





