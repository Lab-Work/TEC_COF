%% Numerical comparison
% This script compares the discrete merge solver and the exact merge
% solver.

% profile on

clearvars

% iteratively run the two solvers at different resolution. For each
% resolution, run 10 tiems and take the average as the mean computation
% time.

% global_num_points = 10:10:5000;   % 100 corresponds to 10 m and 0.3 s
global_num_points = [200];
global_res_counter = 1;

global_exact_timer = zeros(10, length(global_num_points));
global_discrete_timer = zeros(10, length(global_num_points));
global_iter_counter = 1;

for global_pt = 1:length(global_num_points)
    
    for global_iter = 1:1
        
        % clear all parameters to start over
        clearvars -except global_*
        
        % run the two solvers and compute the solution
        %% Configure parameters   
        run_discrete_solver = true;
        run_exact_solver = true;
        
        plot_discrete = false;
        plot_exact = false;
        plot_comparison = true;
        
        
        %% Set the resolution
        % Set the resolution for solving the HJ PDE on each link:
        dx_res_discrete = 1000/global_num_points(global_pt);     % meters
        dt_res_discrete = 30/global_num_points(global_pt);  % seconds
        
        dx_res_exact = 1000/global_num_points(global_pt);       % meters
        dt_res_exact = 30/global_num_points(global_pt);;     % seconds
        
        t_horizon_start = 0;
        sim_steps = 5;
        step_length = 6;       % seconds
        t_horizon_end = t_horizon_start + sim_steps*step_length;
        
        % sol_time = t_horizon_start + dt_res_exact:dt_res_exact:t_horizon_end;  % compute the exact solution at time 30
        
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
        default_para.km_pl = default_para.kc_pl*(default_para.w-default_para.vf)...
            /default_para.w;
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
        % mimic a connection
        % rho_link1 = [1, 1, 1, 1, 1]';
        % rho_link2 = [0, 0, 0, 0, 0]';
        % rho_link3 = [1, 1, 1, 1]';
        
        % with forward and backward shocks
        rho_link1 = [1, 1, 1, 1, 0]';
        rho_link2 = [1, 1, 1, 1, 0]';
        rho_link3 = [4, 1, 1, 1, 1]';
        
        
        %% Define the network boundary conditions
        % Defined as piecewise constant flows on an even time grid.
        % With forward and backward shocks
        q_us_link1 = [1, 0.6, 0.8, 0.8, 0.8]';
        q_us_link2 = [1, 0.6, 0.8, 0.8, 0.8]';
        q_ds_link3 = [0.3, 0.3, 0.3, 1, 1]';
        
        % mimic a connection. Put link 2 all 0
        % q_us_link1 = [1, 0.3, 1, 1, 1]';
        % q_us_link2 = [0, 0, 0, 0, 0]';
        % q_ds_link3 = [1, 1, 1, 1, 1]';
        
        T_BC = ones(5, 1)*step_length;
        
        
        %% Run exact merge solver
        if run_exact_solver == true
            
            exact_timer_start = now;
            % Configure a few parameters for exact solver
            tol_veh = 0.000001;
            tol_tx = [0.0001; 0.0001];
            max_loops = 6;
            
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
            exact_solver.computeExactInternalBoundaryFlows(T_junc, tol_veh, max_loops );
            
            exact_timer_mid = now;
            
            % compute the exact solution at the time t
            M_tx = struct;
            t = 30;
            for link = exact_solver.link_labels'
                
                linkStr = sprintf('link_%d', link);
                % compute the solution at time
                M_tx.(linkStr) = exact_solver.computeSolutionAtTime( t, link, tol_tx );
            end
            
            exact_timer_end = now;
            
%             fprintf('\nExact solver time with %d points:\n',  global_num_points(global_pt))
%             fprintf('-- Solving the junction: %f s\n',...
%                 (exact_timer_mid-exact_timer_start)*86400);
%             fprintf('-- Solving the solution at 30 s: %f s\n',...
%                 (exact_timer_end-exact_timer_mid)*86400);
%             fprintf('-- Total: %f\n', (exact_timer_end-exact_timer_start)*86400);
            
            % save the computation time into the matrix
            global_exact_timer(global_iter, global_pt) = (exact_timer_end-exact_timer_start)*86400;
            
            % visualize the solution in the entire time space domain
            %     for link = exact_solver.link_labels'
            %         % compute solutions on all grid points
            %         exact_solver.computeAllSolutionsOnGrid(dt_res_exact, dx_res_exact);
            %
            %         if plot_exact == true
            %             exact_solver.plotDensityOnLink(link, dt_res_exact, dx_res_exact);
            %             % plot the exact solution for debugging
            %             hold on
            %             linkStr = sprintf('link_%d', link);
            %             scatter(M_tx.(linkStr)(:,1),...
            %                 exact_solver.network_hwy.(linkStr).para_postkm*1000 -...
            %                 M_tx.(linkStr)(:,2),...
            %                 100,...
            %                 'MarkerEdgeColor',[1 0 0],...
            %                 'MarkerFaceColor',[1 0 0],...
            %                 'LineWidth',1.5);
            %
            %             % plot the current grid at the junction
            %             % for upstream
            %             scatter( exact_solver.network_hwy.(linkStr).T_us_cum(2:end,1),....
            %                 ones( length(exact_solver.network_hwy.(linkStr).T_us), 1)*...
            %                 exact_solver.network_hwy.(linkStr).para_postkm*1000,...
            %                 100,...
            %                 'MarkerEdgeColor',[0 0 0],...
            %                 'MarkerFaceColor',[0 0 0],...
            %                 'LineWidth',1.5);
            %             scatter( exact_solver.network_hwy.(linkStr).T_ds_cum(2:end,1),....
            %                 ones( length(exact_solver.network_hwy.(linkStr).T_ds), 1)*0,...
            %                 100,...
            %                 'MarkerEdgeColor',[0 0 0],...
            %                 'MarkerFaceColor',[0 0 0],...
            %                 'LineWidth',1.5);
            %
            %             hold off
            %         end
            %
            %
            %     end
            
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
            discrete_solver.discretizeNet(dt_res_discrete, dx_res_discrete);
            
            % comptue the solution
            discrete_solver.computeSolution(dt_res_discrete, dx_res_discrete);
            
            discrete_timer_end = now;
            
%             fprintf('\nDiscrete solver time with %d points:\n', global_num_points(global_pt))
%             fprintf('-- Total: %f\n', (discrete_timer_end-discrete_timer_start)*86400);
            
            % save the timer 
            global_discrete_timer(global_iter, global_pt) = ...
                (discrete_timer_end-discrete_timer_start)*86400;
            
            if plot_discrete == true
                % Visualizing the solution
                for link = discrete_solver.link_labels'
                    discrete_solver.plotDensityOnLink(link, dt_res_discrete, dx_res_discrete);
                end
            end
            
        end
        
        
        %% plot the comparison of the solution at time t
        if plot_comparison == true
            
            for link = exact_solver.link_labels'
                % plot each link
                linkStr = sprintf('link_%d', link);
                
                scrsz = get(0,'ScreenSize');
                figure('Position',[1 1 scrsz(3) scrsz(4)/2]);
                
                % plot exact solution
                plot( M_tx.(linkStr)(:,2), M_tx.(linkStr)(:,3), 'b', 'Linewidth', 2);
                hold on
                
                % plot approximated solution
                plot(0:dx_res_discrete:discrete_solver.network_hwy.(linkStr).para_postkm*1000,...
                    discrete_solver.M.(linkStr)(:,end), 'r--.', 'Linewidth', 2)
                
                legend('Exact solution', 'Approximated solution');
                title(sprintf('Solution on Link %d at t = %d', link, t), 'fontsize', 30);
                xlabel({'space (m)'},'fontsize',24);
                ylabel({'Solution'},'fontsize',24);
                
                grid on
                hold off
            end
            
            
        end
        
        
        
        
        
        
        
        
        
        
    end
    
    fprintf('\n Solving grid point: %d\n', global_num_points(global_pt))
    
end

% profile viewer
        % p = profile('info');
        % profsave(p,'./result/matlab_profiler_results');

% plot the result
% scrsz = get(0,'ScreenSize');
% figure('Position',[1 1 scrsz(3) scrsz(4)]);
% mean_exact_time = mean(global_exact_timer(:,1:200));
% mean_discrete_time = mean(global_discrete_timer(:,1:200));
% 
% plot(mean_exact_time, 'g', 'Linewidth', 2);
% hold on
% plot(mean_discrete_time, 'r', 'Linewidth', 2);
% grid on
% 
% legend('Exact solutions', 'Approximate solutions', 'fontsize', 24);
% 
% title(linkStr, 'fontsize', 30);
% xlabel({'Number of cells on each link'},'fontsize',24);
% ylabel({'Computation time (s)'},'fontsize',24);


















