%% Numerical comparison
% This script compares the discrete merge solver and the exact merge
% solver.


%% Configure parameters

clearvars -except dbg

run_discrete_solver = true;
run_exact_solver = true;


%%
% Set the resolution for solving the HJ PDE on each link:
dx_res = 1; % meters
dt_res = 0.03; % seconds

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


%% 
% The freeway link 1 & 2 merge to freeway link 3, each with the link length in km:
len_link1 = 1;       % km
len_link2 = 1;
len_link3 = 1;

% number of lanes on each link
num_lanes1 = 2;
num_lanes2 = 1;
num_lanes3 = 2;


%% Define the initial conditions
% Even spacing, normalized to the critical density on each link
rho_link1 = [1, 1, 1, 0, 0]';
rho_link2 = [1, 1, 1, 0, 0]';
rho_link3 = [1, 1, 1, 1, 1]';

%% Define the network boundary conditions
% Defined as piecewise constant flows on an even time grid.
% slower shockwaves inside the link
q_us_link1 = [1, 0.6, 0.8, 0.8, 0.8]';
q_us_link2 = [1, 0.6, 0.8, 0.8, 0.8]';

% more intersection at the boudnary 
% q_us_link1 = [1, 0, 1, 1, 1]';
% q_us_link2 = [1, 0, 1, 1, 1]';
q_ds_link3 = [1, 0, 1, 0, 1]';
T_BC = ones(5, 1)*step_length;


%% Run exact merge solver
if run_exact_solver == true
    
    exact_timer_start = now;
    % Configure a few parameters for exact solver
    tol_veh = 0.1;
    tol_tx = [0.01; 0.1];
    sol_time = 30;  % compute the exact solution at time 30
    
    exact_solver = exactMergeSolver(t_horizon_start, t_horizon_end);
    
    % build the network
    exact_solver.addLink(1, default_para, num_lanes1, len_link1, 'freeway');
    exact_solver.addLink(2, default_para, num_lanes2, len_link2, 'freeway');
    exact_solver.addLink(3, default_para, num_lanes3, len_link3, 'freeway');
    exact_solver.addJunc( 1, [1;2], 3, 'merge', [2;1] );
    
    % set the initial and boundary conditions
    % convert normalized density to veh/m
    Ini.link_1.IC = exact_solver.network_hwy.link_1.para_kc*rho_link1;
    Ini.link_2.IC = exact_solver.network_hwy.link_2.para_kc*rho_link2;
    Ini.link_3.IC = exact_solver.network_hwy.link_3.para_kc*rho_link3;
    exact_solver.setInitialCon(Ini);
    
    % set the network boundary conditions
    % convert the normalized flow to veh/s
    q1_us_data = exact_solver.network_hwy.link_1.para_qmax*q_us_link1;
    q2_us_data = exact_solver.network_hwy.link_2.para_qmax*q_us_link2;
    q3_ds_data = exact_solver.network_hwy.link_3.para_qmax*q_ds_link3;
    
    exact_solver.setBoundaryConForLink(1, q1_us_data, [], T_BC, []);
    exact_solver.setBoundaryConForLink(2, q2_us_data, [], T_BC, []);
    exact_solver.setBoundaryConForLink(3, [], q3_ds_data, [], T_BC);
    
    % compute the exact internal boundary flow solution
    T_junc = T_BC;  % initial grid 
    exact_solver.computeExactInternalBoundaryFlows(T_junc, tol_veh );
    
    exact_timer_mid = now;
    
    % compute solutions on all grid points
    exact_solver.computeAllSolutionsOnGrid(dt_res, dx_res);
    
    exact_timer_end = now;
    
    % visualize the solution in the entire time space domain
    for link = exact_solver.link_labels'
        exact_solver.plotDensityOnLink(link, dt_res, dx_res);
        
        % compute the solution at time
        M_tx = exact_solver.computeSolutionAtTime( sol_time, link, tol_tx );
        
        % plot the exact solution for debugging
        hold on
        linkStr = sprintf('link_%d', link);
        scatter(M_tx(:,1),...
            exact_solver.network_hwy.(linkStr).para_postkm*1000 -...
                M_tx(:,2),...
            100,...
            'MarkerEdgeColor',[1 1 1],...
            'MarkerFaceColor',[1 1 1],...
            'LineWidth',1.5);
        hold off
    end
    
    fprintf('\nExact solver time:\n')
    fprintf('-- Solving the junction: %f s\n',...
                (exact_timer_mid-exact_timer_start)*86400);
    fprintf('-- Solving the states: %f s\n',...
                (exact_timer_end-exact_timer_mid)*86400);
    fprintf('-- Total: %f\n', (exact_timer_end-exact_timer_start)*86400);

    
end
 


%% Run discrete merge solver
if run_discrete_solver == true
    
    discrete_timer_start = now;
    
    discrete_solver = discreteMergeSolver(t_horizon_start, t_horizon_end);
    
    % build the network
    discrete_solver.addLink(1, default_para, num_lanes1, len_link1, 'freeway');
    discrete_solver.addLink(2, default_para, num_lanes2, len_link2, 'freeway');
    discrete_solver.addLink(3, default_para, num_lanes3, len_link3, 'freeway');
    discrete_solver.addJunc( 1, [1;2], 3, 'merge', [2;1] );
    
    % set the initial and boundary conditions
    % convert normalized density to veh/m
    Ini.link_1.IC = discrete_solver.network_hwy.link_1.para_kc*rho_link1;
    Ini.link_2.IC = discrete_solver.network_hwy.link_2.para_kc*rho_link2;
    Ini.link_3.IC = discrete_solver.network_hwy.link_3.para_kc*rho_link3;
    discrete_solver.setInitialCon(Ini);
    
    % set the network boundary conditions
    % convert the normalized flow to veh/s
    q1_us_data = discrete_solver.network_hwy.link_1.para_qmax*q_us_link1;
    q2_us_data = discrete_solver.network_hwy.link_2.para_qmax*q_us_link2;
    q3_ds_data = discrete_solver.network_hwy.link_3.para_qmax*q_ds_link3;
    
    discrete_solver.setBoundaryConForLink(1, q1_us_data, [], T_BC, []);
    discrete_solver.setBoundaryConForLink(2, q2_us_data, [], T_BC, []);
    discrete_solver.setBoundaryConForLink(3, [], q3_ds_data, [], T_BC);
    
    % discretize the network into cells and steps
    discrete_solver.discretizeNet(dt_res, dx_res);
    
    % comptue the solution
    discrete_solver.computeSolution(dt_res, dx_res);
    
    discrete_timer_end = now;
    
    fprintf('\nDiscrete solver time:\n')
    fprintf('-- Total: %f\n', (discrete_timer_end-discrete_timer_start)*86400);
    
    % Visualizing the solution
    for link = discrete_solver.link_labels'
        discrete_solver.plotDensityOnLink(link, dt_res, dx_res);
    end
    
end




























