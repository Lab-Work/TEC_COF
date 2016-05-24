% This example computes the solution to first order conservation law on 
% networks using Godunove scheme. 
% The state is now the vehicle ID.
% We refer to the following paper for the scheme.
% " Intersection modeling using a convergent scheme based on
% Hamilton-Jacobi Equation" Lebacque, 2012.


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
dx_res = 1; % meters
dt_res = 1; % seconds

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
solver = discreteScheme;
solver.addLink(1, default_para, 1, len_link1, 'freeway');
solver.addLink(2, default_para, 1, len_link2, 'freeway');
% 1.5 lanes simulate a limited capacity while still using the default_para
solver.addLink(3, default_para, 1.5, len_link3, 'freeway');   

T_init_grid = ones(sim_steps,1)*step_length;
solver.addJunc(1, [1, 2]', 3, 'merge', [1; 1], T_init_grid);

%% Set the initial and boundary conditions
% Set initial conditoins randomly or staticly.
% Initial traffic density is initialized constants with even
% discretization.

Ini.link_1.IC = solver.network_hwy.link_1.para_kc*[1, 1, 1, 0, 0]';
Ini.link_2.IC = solver.network_hwy.link_2.para_kc*[1, 1, 1, 0, 0]';
Ini.link_3.IC = solver.network_hwy.link_3.para_kc*[1, 1, 1, 1, 1]';

solver.setInitialCon(Ini);

%%
% Set the boundary condition. The downstream boundary flow is not set to
% aviod infeasibility.
q1_us_data = solver.network_hwy.link_1.para_qmax*[1, 0, 1, 1, 1]';
q2_us_data = solver.network_hwy.link_2.para_qmax*[1, 0, 1, 1, 1]';
q3_ds_data = solver.network_hwy.link_3.para_qmax*[1, 0, 1, 0, 1]';


%% Discretize the network into cells and steps
solver.discretizeNet(dt_res, dx_res);


%% Compute the solution
solver.computeSolution();


%% Plot the solution
solver.plotSolution();





















