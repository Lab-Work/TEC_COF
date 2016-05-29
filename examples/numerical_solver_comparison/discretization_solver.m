% This example computes the solution to first order conservation law on 
% networks using Godunove scheme. 
% The state is now the vehicle ID.
% We refer to the following paper for the scheme.
% " Intersection modeling using a convergent scheme based on
% Hamilton-Jacobi Equation" Lebacque, 2012.


%% Configure parameters

clearvars -except dbg

% profile on

timing_start = now;

%%
% Set the resolution for solving the HJ PDE on each link:
dx_res = 50; % meters
dt_res = 1.5; % seconds

t_horizon_start = 0;
sim_steps = 5;
step_length = 30;   % seconds
t_horizon_end = sim_steps*step_length;

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
len_link1 = 1;       % km
len_link2 = 1;
len_link3 = 1;


%% 
% Define the network as an initNetowrk object. Refer to the iniNetwork doc
% for details.
solver = discreteScheme(t_horizon_start, t_horizon_end);
solver.addLink(1, default_para, 2, len_link1, 'freeway');
solver.addLink(2, default_para, 1, len_link2, 'freeway');
% 1.5 lanes simulate a limited capacity while still using the default_para
solver.addLink(3, default_para, 2, len_link3, 'freeway');   

T_init_grid = ones(sim_steps,1)*step_length;
solver.addJunc(1, [1, 2]', 3, 'merge', [2; 1], T_init_grid);

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
q1_us_data = solver.network_hwy.link_1.para_qmax*[1, 0.6, 0.8, 0.8, 0.8]';
q2_us_data = solver.network_hwy.link_2.para_qmax*[1, 0.6, 0.8, 0.8, 0.8]';
q3_ds_data = solver.network_hwy.link_3.para_qmax*[1, 0, 1, 0, 1]';

solver.setBoundaryConForLink(1, q1_us_data, [], T_init_grid, T_init_grid);
solver.setBoundaryConForLink(2, q2_us_data, [], T_init_grid, T_init_grid);
solver.setBoundaryConForLink(3, [], q3_ds_data, T_init_grid, T_init_grid);


%% Discretize the network into cells and steps
solver.discretizeNet(dt_res, dx_res);


%% Compute the solution
solver.computeSolution(dt_res, dx_res);


timing_end = now;

fprintf('Total time: %f\n', (timing_end-timing_start)*86400);



%% Plot the solution
for link = solver.link_labels'
    solver.plotDensityOnLink(link, dt_res, dx_res);
end




















