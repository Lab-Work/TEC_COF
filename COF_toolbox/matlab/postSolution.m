

classdef postSolution < handle
    % This class analyzes and visualizes the solution from CP
    % Yanning Li, Feb 2016
    %   - It integrates the Berkeley HJ PDE solver to solve and visualize
    %     the traffic states on each link.
    %   - It analyzes if the internal boundary flow is admissible, and
    %     updates the boundary grid by locating the wavefront intersection.
    
    properties (Access = protected)
        
        x;  % column vector, solution obtained from the CP
        dv_index;   % save a copy of the index of decision variables from CP
        net;        % save a copy of the network object

        t_start_sim;    % start time of the simulated period, should be 0
        t_end_sim;  % end time of the simulated period
        
        fd; % struct, .(linkStr), fundamental diagram
        dx_res; % meters, space mesh resolution for single link HJ PDE solver
        dt_res; % seconds, time mesh resolution for single link HJ PDE solver
        x_mesh_m; % struct, .(linkStr), space mesh for single link HJ PDE solver
        t_mesh_s; % struct, .(linkStr), time mesh for single link HJ PDE solver
        k;  % struct, .(linkStr), contains the density solution for each link
        N;  % struct, .(linkStr), contains the car ID solution for each link
        
        t_ref;  % the reference point for sampling M.
                %  Note, this is used for determing which boundary conditions 
                %  we should use to compute the sending and receiving function.
        
        hard_queue_limit;   % struct, .(linkStr): in meters from downstream boudnary; just for visualization
        soft_queue_limit;   % struct, .(linkStr): in meters from downstream boudnary; just for visualization

        max_search_depth;   % the maximal seartch depth of the recursive searching the shockwave intersections
        
    end
    
    methods (Access = public)
        
        %===============================================================
        function self=postSolution(x, net, dv_index, end_time, dx_res, dt_res, ...
                        hard_queue_limit, soft_queue_limit)
           % constructs a postSolution object
           % input:
           %       x: float, dv_index_max x 1, solution from cplex
           %       net: the network object constructed by initNetwork
           %       dv_index: the decision varialbe struct
           %       end_time: end time of the simulation
           %       dx_res/dt_res: float, resolution for visualization
           %       hard_queue_limit: struct, .(linkStr) limit of queue in meters
           %       soft_queue_limit: struct, .(linkStr) limit of queue in meters
            self.x = x;
            self.dv_index = dv_index;
            self.net = net;
            self.t_start_sim = 0;
            self.t_end_sim = end_time;
            self.dx_res = dx_res;
            self.dt_res = dt_res;
            self.hard_queue_limit = hard_queue_limit;
            self.soft_queue_limit = soft_queue_limit;
            
            % by default, the maximal search depth is 10
            self.max_search_depth = 10;


        end


        %===============================================================
        function estimateState(self)
            % This function uses Berkeley HJ PDE solver to compute the traffic states on a single link
            % solve each link
            for link = self.net.link_labels'
                
                linkStr = sprintf('link_%d',link);
                
                % extract the computed result 
                q_in = self.x(self.dv_index.(linkStr).upstream(1,1): ...
                    self.dv_index.(linkStr).upstream(2,1));
                q_out = self.x(self.dv_index.(linkStr).downstream(1,1): ...
                    self.dv_index.(linkStr).downstream(2,1));
                p_ini = self.x(self.dv_index.(linkStr).initial(1,1): ...
                    self.dv_index.(linkStr).initial(2,1));
                
                % creat a triangular fundamental diagram
                self.fd.(linkStr) = LH_Tfd(self.net.network_hwy.(linkStr).para_vf,...
                    self.net.network_hwy.(linkStr).para_w,...
                    self.net.network_hwy.(linkStr).para_km);
                
                us_position = 0;
                ds_position = self.net.network_hwy.(linkStr).para_postkm*1000;
                pbEnv = LH_general(self.fd.(linkStr),us_position,ds_position);
                
                
                %===========================================
                % extract initial segment vector
                ini_seg = self.net.network_hwy.(linkStr).X_grid_cum;
                
                % Berkeley toolbox need ini_seg and p_ini to be a row
                % vector
                if ~isrow(ini_seg)
                    ini_seg = ini_seg';
                end
                if ~isrow(p_ini)
                    p_ini = p_ini';
                end
                
                pbEnv.setIniDens(ini_seg,p_ini);
                
                % set upstream boundary condition
                time_grid = self.net.network_hwy.(linkStr).T_us_cum;
                if ~isrow(time_grid)
                    time_grid = time_grid';
                end
                if ~isrow(q_in)
                    q_in = q_in';
                end
                pbEnv.setUsFlows(time_grid,q_in);
                
                % set downstream boundary condition
                time_grid = self.net.network_hwy.(linkStr).T_ds_cum;
                if ~isrow(time_grid)
                    time_grid = time_grid';
                end
                if ~isrow(q_out)
                    q_out = q_out';
                end
                pbEnv.setDsFlows(time_grid,q_out);
                
                %===========================================
                % specify resolution
                nx = floor((ds_position)/self.dx_res);
                dx = (ds_position)/nx;
                
                self.x_mesh_m.(linkStr) = 0:dx:ds_position;
                self.t_mesh_s = self.t_start_sim:self.dt_res:self.t_end_sim;
                
                xValues = ones(size(self.t_mesh_s'))*(self.x_mesh_m.(linkStr));
                tValues = self.t_mesh_s' * ones(size(self.x_mesh_m.(linkStr)));
                
                %===========================================
                % compute Moskowitz
                result = pbEnv.explSol(tValues,xValues);
                
                self.N.(linkStr) = result{1};
                activeComp = result{2};
                self.k.(linkStr) = pbEnv.density(tValues,xValues,activeComp);
            end
            
        end
        
        
        %===============================================================
        function plotLinks(self, links)
            % Plot the density estimate for links. estimateState() function must be called first.
            % input:
            %       links: a column vector of link labels that we want to plot
            
            if strcmp(links,'all')
                links = self.net.link_labels;
            end
            
            
            for link = links'
                
                linkStr = sprintf('link_%d',link);
                T_us_cum = self.net.network_hwy.(linkStr).T_us_cum;
                T_ds_cum = self.net.network_hwy.(linkStr).T_ds_cum;
                len_link = self.net.network_hwy.(linkStr).para_postkm*1000;
                
                %===========================================
                %Transformation for better color presentation
                %kc=>0.5km, km=>km
                k_c_tmp = self.net.network_hwy.(linkStr).para_kc;
                k_m_tmp = self.net.network_hwy.(linkStr).para_km;
                k_trans = self.mapping(self.k.(linkStr), [0 k_c_tmp; k_c_tmp k_m_tmp],...
                    [0 0.5*k_m_tmp; 0.5*k_m_tmp k_m_tmp]);
                
                %===========================================
                scrsz = get(0,'ScreenSize');
                figure('Position',[1 1 scrsz(3) scrsz(4)]);
                title(sprintf('Link No. %d',link),'fontsize',24);
                
                colormap jet
                
                [~, ~] = LH_plot2D(self.t_mesh_s, self.x_mesh_m.(linkStr),...
                    self.N.(linkStr),k_trans, self.fd.(linkStr));
                
                hold on
                % Plot the discretization markers at the upstream
                for i = 1:length(T_us_cum)
                    plot([T_us_cum(i) T_us_cum(i)],...
                         [0 len_link/50],...
                        'k','LineWidth',2);
                end
                
                % plot the discretization markers at the downstream
                for i = 1:length(T_ds_cum)
                    plot([T_ds_cum(i) T_ds_cum(i)],...
                        [len_link-len_link/50 len_link],...
                        'k','LineWidth',2);
                end
                hold off

                set(gca,'fontsize',20)
                xlabel({'time (s)'},'fontsize',24);
                ylabel({'space (m)'},'fontsize',24);
                
            end
               
        end
        
        
        %===============================================================
        function plotJuncs(self, juncs, title_str)
            % plot the density estimates for the links connected to this junction
            %   - connection: plots an density map with two links concatenated
            %   - merge: plots a 1x2 subplots each concatenates two links: 1->3, 2->3
            %   - diverge: plots a 1x2 subplots each concatenates two links: 1->2, 1->3
            % input:
            %       juncs: the column vector of junction labels we want to
            %           plot, each junction will be plot in one seperate figure
            %           or 'all' plots all junctions
            %       title_str: string with some info to put in the title
            
            if strcmp(juncs,'all')
                juncs = self.net.junc_labels;
            end
            
            for junc = juncs'
                
                juncStr = sprintf('junc_%d',junc);
                
                if strcmp(self.net.network_junc.(juncStr).type_junc,'merge') ||...
                   strcmp(self.net.network_junc.(juncStr).type_junc,'onrampjunc') 
                    
                    %===========================================
                    % compute the total through flow and penalty for
                    % not following the priority rule
                    linkStr = sprintf('link_%d', self.net.network_junc.(juncStr).inlabel(1));
                    q_s1 = self.x(self.dv_index.(linkStr).downstream(1,1):...
                                 self.dv_index.(linkStr).downstream(2,1));
                             
                    linkStr = sprintf('link_%d', self.net.network_junc.(juncStr).inlabel(2));
                    q_s2 = self.x(self.dv_index.(linkStr).downstream(1,1):...
                                 self.dv_index.(linkStr).downstream(2,1)); 
                    
                    % the total throughput flow at the junction
                    tt_flow = sum((q_s1+q_s2).*self.net.network_junc.(juncStr).T);
                    
                    % L1 norm penalty for not following the priority rule
                    tt_pen = sum( abs(q_s1 - q_s2*self.net.network_junc.(juncStr).ratio(1)/...
                        self.net.network_junc.(juncStr).ratio(2)).*self.net.network_junc.(juncStr).T );
                    
                    %===========================================
                    % Normalize density
                    % kc=>0.5, km=>1
                    link1 = sprintf('link_%d',self.net.network_junc.(juncStr).inlabel(1));
                    k_c_tmp = self.net.network_hwy.(link1).para_kc;
                    k_m_tmp = self.net.network_hwy.(link1).para_km;
                    kNorm1 = self.mapping(self.k.(link1),...
                        [0 k_c_tmp+0.00001; k_c_tmp+0.00001 k_m_tmp],...
                        [0 0.5; 0.5 1]);
                    
                    link2 = sprintf('link_%d',self.net.network_junc.(juncStr).inlabel(2));
                    k_c_tmp = self.net.network_hwy.(link2).para_kc;
                    k_m_tmp = self.net.network_hwy.(link2).para_km;
                    kNorm2 = self.mapping(self.k.(link2),...
                        [0 k_c_tmp+0.00001; k_c_tmp+0.00001 k_m_tmp],...
                        [0 0.5; 0.5 1]);
                    
                    link3 = sprintf('link_%d',self.net.network_junc.(juncStr).outlabel);
                    k_c_tmp = self.net.network_hwy.(link3).para_kc;
                    k_m_tmp = self.net.network_hwy.(link3).para_km;
                    kNorm3 = self.mapping(self.k.(link3),...
                        [0 k_c_tmp+0.00001; k_c_tmp+0.00001 k_m_tmp],...
                        [0 0.5; 0.5 1]);
                    
                    kLeft = [kNorm1 kNorm3];
                    kRight = [kNorm2 kNorm3];
                    
                    NLeft = [self.N.(link1) self.N.(link3)];
                    NRight = [self.N.(link2) self.N.(link3)];
                    
                    xScaleLeft = [self.x_mesh_m.(link1)...
                        self.x_mesh_m.(link3) + max(self.x_mesh_m.(link1)) + self.dx_res/1000 ];
                    xScaleRight = [self.x_mesh_m.(link2)...
                        self.x_mesh_m.(link3) + max(self.x_mesh_m.(link2)) + self.dx_res/1000];
                    
                    %===========================================
                    scrsz = get(0,'ScreenSize');
                    figure('Position',[1 1 scrsz(3) scrsz(4)]);
                    %title(sprintf('Link No. %d',link(i)),'fontsize',24);
                    
                    colormap jet
                    
                    [ax1, ax2] = LH_plot2D_junc([self.net.network_hwy.(link1).para_postkm*1000,...
                        self.net.network_hwy.(link2).para_postkm*1000],...
                        self.net.network_junc.(juncStr).T_cum',...
                        length(self.net.network_junc.(juncStr).T), ...
                        self.t_mesh_s, xScaleLeft, xScaleRight,...
                        NLeft, kLeft, NRight, kRight, self.net.network_junc.(juncStr));
                    h = suptitle(sprintf('%s\n Total Flow: %f;\n Number of steps: %d; Penalty: %f',...
                        title_str, tt_flow, length(self.net.network_junc.(juncStr).T), tt_pen));
                    set(h,'FontSize',20)
                    
                    % plot additional hard queue limit solid lines and the
                    % soft queue limit in dash lines
                    t_line_hard = [self.net.network_junc.(juncStr).T_cum(1),...
                                  self.net.network_junc.(juncStr).T_cum(end)]';
                    t_line_soft = self.t_start_sim:30:self.t_end_sim;
                    if t_line_soft(end) ~= self.t_end_sim
                        t_line_soft = [t_line_soft'; self.t_end_sim];
                    else
                        t_line_soft = t_line_soft';
                    end
                              
                    hold(ax1, 'on');
                    if isfield(self.hard_queue_limit, link1)
                        x_queue = self.net.network_hwy.(link1).para_postkm*1000 - self.hard_queue_limit.(link1);
                        plot(ax1, t_line_hard, [x_queue, x_queue]', 'r', 'LineWidth', 3)
                    end
                    if isfield(self.hard_queue_limit, link3)
                        x_queue = self.net.network_hwy.(link1).para_postkm*1000 + ...
                                  self.net.network_hwy.(link3).para_postkm*1000 - self.hard_queue_limit.(link3);
                        plot(ax1, t_line_hard, [x_queue, x_queue]', 'r', 'LineWidth', 3)
                    end
                    if isfield(self.soft_queue_limit, link1)
                        x_queue = self.net.network_hwy.(link1).para_postkm*1000 - self.soft_queue_limit.(link1);
                        plot(ax1, t_line_soft, x_queue*ones(length(t_line_soft), 1), '+r--', 'LineWidth', 3)
                    end
                    if isfield(self.soft_queue_limit, link3)
                        x_queue = self.net.network_hwy.(link1).para_postkm*1000 + ...
                                  self.net.network_hwy.(link3).para_postkm*1000 - self.soft_queue_limit.(link3);
                        plot(ax1,  t_line_soft, x_queue*ones(length(t_line_soft), 1), '*r--', 'LineWidth', 3)
                    end
                    hold(ax1, 'off');
                    
                    hold(ax2,'on');
                    if isfield(self.hard_queue_limit, link2)
                        x_queue = self.net.network_hwy.(link2).para_postkm*1000 - self.hard_queue_limit.(link2);
                        plot(ax2, t_line_hard, [x_queue, x_queue]', 'r', 'LineWidth', 3)
                    end
                    if isfield(self.hard_queue_limit, link3)
                        x_queue = self.net.network_hwy.(link2).para_postkm*1000 + ...
                                  self.net.network_hwy.(link3).para_postkm*1000 - self.hard_queue_limit.(link3);
                        plot(ax2, t_line_hard, [x_queue, x_queue]', 'r', 'LineWidth', 3)
                    end
                    if isfield(self.soft_queue_limit, link2)
                        x_queue = self.net.network_hwy.(link2).para_postkm*1000 - self.soft_queue_limit.(link2);
                        plot(ax2, t_line_soft, x_queue*ones(length(t_line_soft), 1), '+r--', 'LineWidth', 3)
                    end
                    if isfield(self.soft_queue_limit, link3)
                        x_queue = self.net.network_hwy.(link2).para_postkm*1000 + ...
                                  self.net.network_hwy.(link3).para_postkm*1000 - self.soft_queue_limit.(link3);
                        plot(ax2, t_line_soft, x_queue*ones(length(t_line_soft), 1), '*r--', 'LineWidth', 3)
                    end
                    hold(ax2, 'off');
                    
                    
                    
                elseif strcmp(self.net.network_junc.(juncStr).type_junc,'diverge') ||...
                       strcmp(self.net.network_junc.(juncStr).type_junc,'offrampjunc')
                    %===========================================
                    % compute the total through flow and penalty for
                    % not following the distribution rule
                    linkStr = sprintf('link_%d', self.net.network_junc.(juncStr).outlabel(1));
                    q_r1 = self.x(self.dv_index.(linkStr).upstream(1,1):...
                                 self.dv_index.(linkStr).upstream(2,1));
                             
                    linkStr = sprintf('link_%d', self.net.network_junc.(juncStr).outlabel(2));
                    q_r2 = self.x(self.dv_index.(linkStr).upstream(1,1):...
                                 self.dv_index.(linkStr).upstream(2,1)); 
                                       
                    tt_flow = sum((q_r1+q_r2).*self.net.network_junc.(juncStr).T);
                    % L1 norm penalty
                    tt_pen = sum( abs(q_r1 - q_r2*self.net.network_junc.(juncStr).ratio(1)/...
                        self.net.network_junc.(juncStr).ratio(2)).*self.net.network_junc.(juncStr).T );
                    
                    
                    %===========================================
                    %Normalize density
                    %kc=>0.5, km=>1
                    link1 = sprintf('link_%d',self.net.network_junc.(juncStr).inlabel);
                    k_c_tmp = self.net.network_hwy.(link1).para_kc;
                    k_m_tmp = self.net.network_hwy.(link1).para_km;
                    kNorm1 = self.mapping(self.k.(link1),...
                        [0 k_c_tmp+0.00001; k_c_tmp+0.00001 k_m_tmp],...
                        [0 0.5; 0.5 1]);
                    
                    link2 = sprintf('link_%d',self.net.network_junc.(juncStr).outlabel(1));
                    k_c_tmp = self.net.network_hwy.(link2).para_kc;
                    k_m_tmp = self.net.network_hwy.(link2).para_km;
                    kNorm2 = self.mapping(self.k.(link2),...
                        [0 k_c_tmp+0.00001; k_c_tmp+0.00001 k_m_tmp],...
                        [0 0.5; 0.5 1]);
                    
                    link3 = sprintf('link_%d',self.net.network_junc.(juncStr).outlabel(2));
                    k_c_tmp = self.net.network_hwy.(link3).para_kc;
                    k_m_tmp = self.net.network_hwy.(link3).para_km;
                    kNorm3 = self.mapping(self.k.(link3),...
                        [0 k_c_tmp+0.00001; k_c_tmp+0.00001 k_m_tmp],...
                        [0 0.5; 0.6 1]);
                    
                    kLeft = [kNorm1 kNorm2];
                    kRight = [kNorm1 kNorm3];
                    
                    NLeft = [self.N.(link1) self.N.(link2)];
                    NRight = [self.N.(link1) self.N.(link3)];
                    
                    xScaleLeft = [self.x_mesh_m.(link1), ...
                        self.x_mesh_m.(link2) + max(self.x_mesh_m.(link1)) + self.dx_res/1000 ];
                    xScaleRight = [self.x_mesh_m.(link1), ...
                        self.x_mesh_m.(link3) + max(self.x_mesh_m.(link1)) + self.dx_res/1000 ];
                    
                    %===========================================
                    scrsz = get(0,'ScreenSize');
                    figure('Position',[1 1 scrsz(3) scrsz(4)]);
                    %title(sprintf('Link No. %d',link(i)),'fontsize',24);
                    
                    colormap jet
                    
                    [~, ~] = LH_plot2D_junc([self.net.network_hwy.(link1).para_postkm*1000,...
                        self.net.network_hwy.(link1).para_postkm*1000],...
                        self.net.network_junc.(juncStr).T_cum',...
                        length(self.net.network_junc.(juncStr).T), ...
                        self.t_mesh_s, xScaleLeft, xScaleRight,...
                        NLeft, kLeft, NRight, kRight,self.net.network_junc.(juncStr));
                    h = suptitle(sprintf('%s\n Total Flow: %f;\n Number of steps: %d; Penalty: %f\n',...
                        title_str, tt_flow, length(self.net.network_junc.(juncStr).T)), tt_pen);
                    set(h,'FontSize',24)
                    
                elseif strcmp(self.net.network_junc.(juncStr).type_junc,'connection')
                    %===========================================
                    linkStr = sprintf('link_%d', self.net.network_junc.(juncStr).outlabel);
                    q_thru = self.x( self.dv_index.(linkStr).upstream(1,1):...
                                    self.dv_index.(linkStr).upstream(2,1));
                    
                    tt_flow = sum((q_thru).*self.net.network_junc.(juncStr).T);
                    
                    %===========================================
                    %Normalize density
                    %kc=>0.5, km=>1
                    link1 = sprintf('link_%d',self.net.network_junc.(juncStr).inlabel);
                    k_c_tmp = self.net.network_hwy.(link1).para_kc;
                    k_m_tmp = self.net.network_hwy.(link1).para_km;
                    kNorm1 = self.mapping(self.k.(link1),...
                        [0 k_c_tmp+0.00001; k_c_tmp+0.00001 k_m_tmp],...
                        [0 0.5; 0.5 1]);
                    
                    link2 = sprintf('link_%d',self.net.network_junc.(juncStr).outlabel);
                    k_c_tmp = self.net.network_hwy.(link2).para_kc;
                    k_m_tmp = self.net.network_hwy.(link2).para_km;
                    kNorm2 = self.mapping(self.k.(link2),...
                        [0 k_c_tmp+0.00001; k_c_tmp+0.00001 k_m_tmp],...
                        [0 0.5; 0.5 1]);
                    
                    kComb = [kNorm1 kNorm2];
                    
                    NComb = [self.N.(link1) self.N.(link2)];
                    
                    xScaleComb = [self.x_mesh_m.(link1)...
                        self.x_mesh_m.(link2) + max(self.x_mesh_m.(link1)) + self.dx_res/1000];
                    
                    %===========================================
                    scrsz = get(0,'ScreenSize');
                    figure('Position',[1 1 scrsz(3) scrsz(4)]);
                    %title(sprintf('Link No. %d',link(i)),'fontsize',24);
                    
                    colormap jet
                    
                    [~, ~] = LH_plot2D_junc(self.net.network_hwy.(link1).para_postkm*1000,...
                        self.net.network_junc.(juncStr).T_cum',...
                        length(self.net.network_junc.(juncStr).T), ... 
                        self.t_mesh_s, xScaleComb,...
                        NComb, kComb,self.net.network_junc.(juncStr));
                    h = suptitle(sprintf('Total Flow: %f\n Number of steps: %d\n',...
                        tt_flow, length(self.net.network_junc.(juncStr).T)));
                    set(h,'FontSize',24)
                end
                
            end
            
        end     
        
        
        %===============================================================
        function [TF, steps] = checkEntropy(self, entropyTol)
            % Check if the solution is admissible.
            %   - uses entropicSolutionAtJunc to compute the admissible
            %       solution
            %   - uses sending and receiving functions to check if the flow
            %       is maximized while minimizing the priority ratio.
            % input: 
            %       net: the network object 
            %       entropy_tol: cars, the tolerance for entropy, for instance, 
            %           if the error to entropic solution is only 1 car, we
            %           consider this solution as an entropic solution and stop
            %           updating the discretization grid.
            % output: 
            %       TF: true or false, this solution is an entropic solution?
            %       steps: struct, .(juncStr), the id of the first non entropic
            %           step
                        
            TF = true;
            steps = struct;
            
            % check each junction
            for junc = self.net.junc_labels'
                
                juncStr = sprintf('junc_%d',junc);
                
                % initialize as empty
                steps.(juncStr) = [];
                
                % the boudnary discretization grid at this junction;
                % they are also saved in the links
                T_grid = self.net.network_junc.(juncStr).T;
                T_cum_grid = self.net.network_junc.(juncStr).T_cum;
                
                % Here consider connection case
                if strcmp(self.net.network_junc.(juncStr).type_junc,'connection')
                    
                    outlink = self.net.network_junc.(juncStr).outlabel;
                    linkStr = sprintf('link_%d', outlink);
                    
                    % extract the thru flow value which is downstream
                    % inflow in the connection case
                    q_thru = self.x(self.dv_index.(linkStr).upstream(1,1):...
                                   self.dv_index.(linkStr).upstream(2,1));
                    
                    % check each step
                    % break once found the first non-entropic step
                    for i = 1: length(T_grid)
                        
                        % set the reference point of this step
                        self.t_ref = T_cum_grid(i);   
                        t_start = T_cum_grid(i);
                        t_end = T_cum_grid(i+1); % end time of this point
                        
                        % Sample the start and end point of this step;
                        % compute the min of sending and receiving function
                        % which should be entropic solution
                        d_M = self.entropicSolutionAtJunc(t_end, junc);
                        
                        if abs(q_thru(i)*T_grid(i)-d_M) > entropyTol
                            % difference greater than the tolerance
                            % This may due to the discretization or we
                            % simply did not add entropic component in the
                            % objective function
                            d_M_SR = [self.sendingFuncAtLink(t_end, inlink);
                                      self.receivingFuncAtLink(t_end, outlink)];
                            t_C = (t_start+t_end)/2;
                            d_M_C = [self.sendingFuncAtLink(t_C, inlink);
                                     self.receivingFuncAtLink(t_C, outlink)];
                            
                            if self.onStraightLine([t_start, t_C, t_end]',[0, d_M_C(1), d_M_SR(1)]') && ...
                               self.onStraightLine([t_start, t_C, t_end]',[0, d_M_C(2), d_M_SR(2)]')     
                                % meaning caused by not setting entropic
                                % solution or traffic control
                                warning('WARNING: Step %d is not entropic\n', i);
                                continue
                            else
                                TF = false;
                                steps.(juncStr) = i;
                                break
                            end
                            
                        else
                            % If the first half sends more and second half
                            % sends less, then the average may be same as
                            % the entropic solution. In this case, the
                            % affected domain of the discretization error
                            % is small. we treat this as entropic
                            continue
                        end
                    end     %end each step
                
                % Here consider onrampjunc case
                elseif strcmp(self.net.network_junc.(juncStr).type_junc,'onrampjunc')
                    
                    % find the upstream freeway, onramp, and the downstream
                    % freeway
                    inlinks = self.net.network_junc.(juncStr).inlabel;
                    outlink = self.net.network_junc.(juncStr).outlabel;
                    linkStr1 = sprintf('link_%d',inlinks(1));
                    linkStr2 = sprintf('link_%d',inlinks(2));
                    
                    if strcmp(self.net.network_hwy.(linkStr1).para_linktype, 'freeway')
                        upFreewayStr = linkStr1;
                        up_link = inlinks(1);
                    else
                        upFreewayStr = linkStr2;
                        up_link = inlinks(2);
                    end
                    
                    % extract the upstream link flow
                    q_up = self.x(self.dv_index.(upFreewayStr).downstream(1,1):...
                                   self.dv_index.(upFreewayStr).downstream(2,1));
                    
                    % check each step
                    % break once found the first non-entropic step
                    for i = 1: length(T_grid)
                        
                        % set the reference point of this step
                        self.t_ref = T_cum_grid(i);   
                        t_start = T_cum_grid(i);
                        t_end = T_cum_grid(i+1); % end time of this point
                        
                        % find out the receiving funciton on the downstream
                        d_M = self.entropicSolutionAtJunc(t_end, junc);
                        
                        if abs(q_up(i)*T_grid(i) - d_M) > entropyTol
                            % difference greater than the tolerance
                            % This may due to the discretization or we
                            % simply did not add entropic component in the
                            % objective function
                            d_M_SR = [self.sendingFuncAtLink(t_end, up_link);
                                      self.receivingFuncAtLink(t_end, outlink)];
                            t_C = (t_start+t_end)/2;
                            d_M_C = [self.sendingFuncAtLink(t_C, up_link);
                                     self.receivingFuncAtLink(t_C, outlink)];
                            
                            if self.onStraightLine([t_start, t_C, t_end]',[0, d_M_C(1), d_M_SR(1)]') && ...
                               self.onStraightLine([t_start, t_C, t_end]',[0, d_M_C(2), d_M_SR(2)]')     
                                % meaning caused by not setting entropic
                                % solution or traffic control
                                warning('WARNING: Step %d is not entropic\n', i);
                                continue
                            else
                                TF = false;
                                steps.(juncStr) = i;
                                break
                            end
                            
                        else
                            % If the first half sends more and second half
                            % sends less, then the average may be same as
                            % the entropic solution. In this case, the
                            % affected domain of the discretization error
                            % is small. we treat this as entropic
                            continue
                        end
                    end     %end each step    
                    
                % Here consider offrampjunc case
                elseif strcmp(self.net.network_junc.(juncStr).type_junc,'offrampjunc')
                    
                    % find the the downstream freeway
                    inlink = self.net.network_junc.(juncStr).inlabel;
                    outlinks = self.net.network_junc.(juncStr).outlabel;
                    linkStr1 = sprintf('link_%d',outlinks(1));
                    linkStr2 = sprintf('link_%d',outlinks(2));
                    
                    if strcmp(self.net.network_hwy.(linkStr1).para_linktype, 'freeway')
                        downFreewayStr = linkStr1;
                        down_link = inlinks(1);
                    else
                        downFreewayStr = linkStr2;
                        down_link = inlinks(2);
                    end
                    
                    % extract the downstream link flow
                    q_down = self.x(self.dv_index.(downFreewayStr).upstream(1,1):...
                                   self.dv_index.(downFreewayStr).upstream(2,1));
                    
                    % check each step
                    % break once found the first non-entropic step
                    for i = 1: length(T_grid)
                        
                        % set the reference point of this step
                        self.t_ref = T_cum_grid(i);   
                        t_start = T_cum_grid(i);
                        t_end = T_cum_grid(i+1); % end time of this point
                        
                        % find out the receiving funciton on the downstream
                        d_M = self.entropicSolutionAtJunc(t_end, junc);
                        
                        if abs(q_down(i)*T_grid(i) - d_M) > entropyTol
                            % difference greater than the tolerance
                            % This may due to the discretization or we
                            % simply did not add entropic component in the
                            % objective function
                            d_M_SR = [self.sendingFuncAtLink(t_end, inlink);
                                      self.receivingFuncAtLink(t_end, down_link)];
                            t_C = (t_start+t_end)/2;
                            d_M_C = [self.sendingFuncAtLink(t_C, inlink);
                                     self.receivingFuncAtLink(t_C, down_link)];
                            
                            if self.onStraightLine([t_start, t_C, t_end]',[0, d_M_C(1), d_M_SR(1)]') && ...
                               self.onStraightLine([t_start, t_C, t_end]',[0, d_M_C(2), d_M_SR(2)]')     
                                % meaning caused by not setting entropic
                                % solution or traffic control
                                warning('WARNING: Step %d is not entropic\n', i);
                                continue
                            else
                                TF = false;
                                steps.(juncStr) = i;
                                break
                            end
                            
                        else
                            % If the first half sends more and second half
                            % sends less, then the average may be same as
                            % the entropic solution. In this case, the
                            % affected domain of the discretization error
                            % is small. we treat this as entropic
                            continue
                        end
                    end     %end each step 
                    
                    
                elseif strcmp(self.net.network_junc.(juncStr).type_junc,'merge')
                    
                    % extract the outflow of the two incoming links
                    inlinks = self.net.network_junc.(juncStr).inlabel;
                    outlink = self.net.network_junc.(juncStr).outlabel;
                    linkStr = sprintf('link_%d',inlinks(1));
                    q_s1 = self.x(self.dv_index.(linkStr).downstream(1,1):...
                                 self.dv_index.(linkStr).downstream(2,1) );
                    linkStr = sprintf('link_%d',inlinks(2));
                    q_s2 = self.x(self.dv_index.(linkStr).downstream(1,1):...
                                 self.dv_index.(linkStr).downstream(2,1) );
                    
                    % check each step
                    for i = 1: length(T_grid)
                        
                        self.t_ref = T_cum_grid(i);   % start time of this point
                        t_start = T_cum_grid(i);
                        t_end = T_cum_grid(i+1); % end time of this point
                       
                        % a 1 x 2 vector true unique solution from the two incoming links
                        d_M = self.entropicSolutionAtJunc( t_end, junc);
                        
                        if abs( q_s1(i)*T_grid(i)-d_M(1) ) > entropyTol ||...
                                abs( q_s2(i)*T_grid(i)-d_M(2) ) > entropyTol
                            % does not match the unique solution from the
                            % sampling approach
                            % the sending and receiving functions on the
                            % links, just to check there is not
                            % intersection on any three links
                            d_M_SR = [self.sendingFuncAtLink(t_end, inlinks(1));
                                      self.sendingFuncAtLink(t_end, inlinks(2));
                                      self.receivingFuncAtLink(t_end, outlink)];
                            t_C = (t_start+t_end)/2;
                            d_M_C = [self.sendingFuncAtLink(t_C, inlinks(1));
                                     self.sendingFuncAtLink(t_C, inlinks(2));
                                     self.receivingFuncAtLink(t_C, outlink)];
                            
                            if self.onStraightLine([t_start, t_C, t_end]',[0, d_M_C(1), d_M_SR(1)]') &&...
                                  self.onStraightLine([t_start, t_C, t_end]',[0, d_M_C(2), d_M_SR(2)]') &&...
                                  self.onStraightLine([t_start, t_C, t_end]',[0, d_M_C(3), d_M_SR(3)]')
                                % meaning caused by not setting entropic
                                % solution or traffic control
                                warning('WARNING: Step %d is not entropic\n', i);
                                continue
                            else
                                TF = false;
                                steps.(juncStr) = i;
                                break
                            end
                            
                        else
                            % If the first half sends more and second half
                            % sends less, then the average may be same as
                            % the entropic solution. Hence double check the
                            % center point.
                            t_C = (t_start+t_end)/2;
                            d_M_C = self.entropicSolutionAtJunc( t_C, junc);
                            
                            if abs( q_s1(i)*T_grid(i)/2 - d_M_C(1) ) > entropyTol ||...
                                    abs( q_s2(i)*T_grid(i)/2 - d_M_C(2) ) > entropyTol
                                TF = false;
                                steps.(juncStr) = i;
                                break
                            else
                                % entropy solutution
                                continue
                            end
                        end
                    end     % end each step
                    
                elseif strcmp(self.net.network_junc.(juncStr).type_junc,'diverge')
                                        
                    % extract the outflow of the two incoming links
                    outlinks = self.net.network_junc.(juncStr).outlabel;
                    inlink = self.net.network_junc.(juncStr).inlabel;
                    linkStr = sprintf('link_%d',outlinks(1));
                    q_r1 = self.x(self.dv_index.(linkStr).upstream(1,1):...
                                 self.dv_index.(linkStr).upstream(2,1) );
                    linkStr = sprintf('link_%d',outlinks(2));
                    q_r2 = self.x(self.dv_index.(linkStr).upstream(1,1):...
                                 self.dv_index.(linkStr).upstream(2,1) );
                    
                    % check each step
                    for i = 1: length(T_grid)
                        
                        self.t_ref = T_cum_grid(i);   % start time of this point
                        t_start = T_cum_grid(i);
                        t_end = T_cum_grid(i+1); % end time of this point
                        
                        % a 1 x 2 vector true unique solution from the two incoming links
                        d_M = self.entropicSolutionAtJunc( t_end, junc);
                        
                        if abs( q_r1(i)*T_grid(i)-d_M(1) ) > entropyTol ||...
                                abs( q_r2(i)*T_grid(i)-d_M(2) ) > entropyTol
                            % does not match the unique solution from the
                            % sampling approach
                            d_M_SR = [self.receivingFuncAtLink(t_end, outlinks(1));
                                      self.receivingFuncAtLink(t_end, outlinks(2));
                                      self.sendingFuncAtLink(t_end, inlink)];
                            t_C = (t_start+t_end)/2;
                            d_M_C = [self.receivingFuncAtLink(t_C, outlinks(1));
                                     self.receivingFuncAtLink(t_C, outlinks(2));
                                     self.sendingFuncAtLink(t_C, inlink);];
                            
                            if self.onStraightLine([t_start, t_C, t_end]',[0, d_M_C(1), d_M_SR(1)]') && ...
                               self.onStraightLine([t_start, t_C, t_end]',[0, d_M_C(2), d_M_SR(2)]') && ...
                               self.onStraightLine([t_start, t_C, t_end]',[0, d_M_C(3), d_M_SR(3)]')
                                % meaning caused by not setting entropic
                                % solution or traffic control
                                warning('WARNING: Step %d is not entropic\n', i);
                                continue
                            else
                                TF = false;
                                steps.(juncStr) = i;
                                break
                            end
                            
                        else
                            % If the first half sends more and second half
                            % sends less, then the average may be same as
                            % the entropic solution. 
                            
                            t_C = (t_start+t_end)/2;
                            d_M_C = self.entropicSolutionAtJunc( t_C, junc);
                            
                            if abs( q_r1(i)*T_grid(i)/2 - d_M_C(1) ) > entropyTol ||...
                                    abs( q_r2(i)*T_grid(i)/2 - d_M_C(2) ) > entropyTol
                                TF = false;
                                steps.(juncStr) = i;
                                break
                            else
                                % entropy solutution
                                continue
                            end
                        end
                    end     % end each step
                    
                end     % end diverge
                
                
            end % end each junction
            
            % steps = unique(steps);
            
        end 
        
        
        %===============================================================
        function T_new_grid = updateTimeDiscretization(self, steps)
            % Updates the the time grid by searching the intersection points
            %   - This function uses the old versions:
            %       searchBoundaryIntersection(), findBoundaryFuncSlope(),
            %   - Find the frontwave intersection point to eliminate the
            %       discretization error.
            % input: 
            %       steps: struct, .(juncStr), the first non-entropic step. 
            %           Note, the update has to be done iteratively, 
            %           earlier steps first
            % output: 
            %       T_new_cum: struct, .(juncStr). The updated time grid. 
                        
            % recursively search uses a binary cut method to locate the
            % intersection point. 
            % Set the following parameters to limit the number of
            % iterations: stop if we think the position is roughtly good.
            searchDepth = 0;
            d_t = 1.0e-3;    % 0.001 second
            
            % the new discretization grid to be returned.
            T_new_grid = struct;
            
            % check each junction
            for junc = self.net.junc_labels'
                
                juncStr = sprintf('junc_%d',junc);
                
                % get the link label at this junction
                inlink = self.net.network_junc.(juncStr).inlabel;
                outlink = self.net.network_junc.(juncStr).outlabel;
                
                T_tmp_cum = self.net.network_junc.(juncStr).T_cum;
                end_time = self.net.network_junc.(juncStr).T_cum(end);
                num_steps = length(self.net.network_junc.(juncStr).T);
                
                % find the start and end time of the next step of the
                % non-entropic step, since the non-entropy could be caused 
                % by the intersection from the following step.
                if ~isempty(steps.(juncStr))
                    
                    nonentropic_step = steps.(juncStr);
                    
                    self.t_ref = self.net.network_junc.(juncStr).T_cum(nonentropic_step);
                    t_left = self.net.network_junc.(juncStr).T_cum(nonentropic_step);
                
                    if nonentropic_step < num_steps
                        t_right = self.net.network_junc.(juncStr).T_cum(nonentropic_step + 2);
                    else
                        t_right = self.net.network_junc.(juncStr).T_cum(nonentropic_step+1);
                    end
                end
                    
                
                if strcmp( self.net.network_junc.(juncStr).type_junc, 'connection')
                    
                    if ~isempty(steps.(juncStr))
                        
                        % search the intersection of shockwaves on incoming links
                        % t_found = searchBoundaryIntersection(obj, t_interval, M_interval, slope, ...
                        %                      link, bound, searchDepth, dt_tol)
                        t_found_s = searchBoundaryIntersection(self, [t_left, t_right]', [0, NaN]', [NaN, NaN]',...
                                            inlink, 'downstream', searchDepth, d_t);
                        % search the intersection of shockwaves on outgoing links
                        t_found_r = searchBoundaryIntersection(self, [t_left, t_right]', [0, NaN]', [NaN, NaN]',...
                                            outlink, 'upstream', searchDepth, d_t);
                        
                        t_found = [t_found_s; t_found_r];
                    else
                        t_found = [];
                    end
                    
                    
                    % Aggregate original step intervals if the flows in two
                    % consecutive steps are not changing
                    % This reduces the number of discretized intervals in
                    % the boudnary without causing discretization error
                    linkStr = sprintf('link_%d', outlink);
                    q_thru = self.x(self.dv_index.(linkStr).upstream(1,1):...
                                   self.dv_index.(linkStr).upstream(2,1));
                    tmp_group = self.groupSameElement(q_thru, inf);
                    T_tmp_cum = [T_tmp_cum(tmp_group(:,1));...
                                 end_time];
                    T_tmp_cum = unique(T_tmp_cum);
                    
                % if the onrampjunc. Then find the boundary intersection
                % point at the upstream freeway and downstream freeway
                elseif strcmp( self.net.network_junc.(juncStr).type_junc, 'onrampjunc')
                    
                    % find the upstream freeway, onramp, and the downstream
                    % freeway
                    inlinks = self.net.network_junc.(juncStr).inlabel;
                    outlink = self.net.network_junc.(juncStr).outlabel;
                    linkStr1 = sprintf('link_%d',inlinks(1));
                    
                    if strcmp(self.net.network_hwy.(linkStr1).para_linktype, 'freeway')
                        up_link = inlinks(1);
                    else
                        up_link = inlinks(2);
                    end
                    
                    if ~isempty(steps.(juncStr))
                        
                        % search the intersection of shockwaves on incoming links
                        % t_found = searchBoundaryIntersection(obj, t_interval, M_interval, slope, ...
                        %                      link, bound, searchDepth, dt_tol)
                        t_found_s = searchBoundaryIntersection(self, [t_left, t_right]', [0, NaN]',...
                                            [NaN, NaN]', up_link, 'downstream', searchDepth, d_t);
                        % search the intersection of shockwaves on outgoing links
                        t_found_r = searchBoundaryIntersection(self, [t_left, t_right]', [0, NaN]',...
                                            [NaN, NaN]', outlink, 'upstream', searchDepth, d_t);
                        
                        t_found = [t_found_s; t_found_r];
                    else
                        t_found = [];
                    end
                    
                    % Aggregate original step intervals if the flows in two
                    % consecutive steps are not changing
                    linkStr = sprintf('link_%d', inlink(1));
                    q_s1 = self.x( self.dv_index.(linkStr).downstream(1,1):...
                                  self.dv_index.(linkStr).downstream(2,1) );
                    linkStr = sprintf('link_%d', inlink(2));
                    q_s2 = self.x( self.dv_index.(linkStr).downstream(1,1):...
                                  self.dv_index.(linkStr).downstream(2,1) );
                    
                    T_tmp_cum_s1 = T_tmp_cum;
                    tmp_group = self.groupSameElement(q_s1, inf);
                    T_tmp_cum_s1 = [ T_tmp_cum_s1(tmp_group(:,1));...
                                     end_time];
                    T_tmp_cum_s2 = T_tmp_cum;
                    tmp_group = self.groupSameElement(q_s2, inf);
                    T_tmp_cum_s2 = [ T_tmp_cum_s2(tmp_group(:,1));...
                                     end_time];
                    
                    T_tmp_cum = unique([T_tmp_cum_s1; T_tmp_cum_s2]);
                    
                % if the offrampjunc. Then find the boundary intersection
                % point at the upstream freeway and downstream freeway
                elseif strcmp( self.net.network_junc.(juncStr).type_junc, 'offrampjunc')
                    
                    % find the upstream freeway, onramp, and the downstream
                    % freeway
                    inlink = self.net.network_junc.(juncStr).inlabel;
                    outlinks = self.net.network_junc.(juncStr).outlabel;
                    linkStr1 = sprintf('link_%d',outlinks(1));
                    
                    if strcmp(self.net.network_hwy.(linkStr1).para_linktype, 'freeway')
                        down_link = outlinks(1);
                    else
                        down_link = outlinks(2);
                    end
                    
                    if ~isempty(steps.(juncStr))
                        
                        % search the intersection of shockwaves on incoming links
                        % t_found = searchBoundaryIntersection(obj, t_interval, M_interval, slope, ...
                        %                      link, bound, searchDepth, dt_tol)
                        t_found_s = searchBoundaryIntersection(self, [t_left, t_right]', [0, NaN]', [NaN, NaN]',...
                                            inlink, 'downstream', searchDepth, d_t);
                        % search the intersection of shockwaves on outgoing links
                        t_found_r = searchBoundaryIntersection(self, [t_left, t_right]', [0, NaN]', [NaN, NaN]',...
                                            down_link, 'upstream', searchDepth, d_t);
                        
                        t_found = [t_found_s; t_found_r];
                    else
                        t_found = [];
                    end
                    
                    % Aggregate original step intervals if the flows in two
                    % consecutive steps are not changing
                    linkStr = sprintf('link_%d', outlink(1));
                    q_r1 = self.x( self.dv_index.(linkStr).upstream(1,1):...
                                  self.dv_index.(linkStr).upstream(2,1) );
                    linkStr = sprintf('link_%d', outlink(2));
                    q_r2 = self.x( self.dv_index.(linkStr).upstream(1,1):...
                                  self.dv_index.(linkStr).upstream(2,1) );
                    
                    T_tmp_cum_r1 = T_tmp_cum;
                    tmp_group = self.groupSameElement(q_r1, inf);
                    T_tmp_cum_r1 = [ T_tmp_cum_r1(tmp_group(:,1));...
                                     end_time];
                    T_tmp_cum_r2 = T_tmp_cum;
                    tmp_group = self.groupSameElement(q_r2, inf);
                    T_tmp_cum_r2 = [ T_tmp_cum_r2(tmp_group(:,1));...
                                     end_time];
                    
                    T_tmp_cum = unique([T_tmp_cum_r1; T_tmp_cum_r2]);
                    
                
                % if merge junction, then find the intersection point 
                % bewteen shock waves and the boundary on three links
                elseif strcmp( self.net.network_junc.(juncStr).type_junc, 'merge')
                    
                    if ~isempty(steps.(juncStr))
                    
                        % corner points found from sending link 1, 2 and
                        % receiving link 3
                        t_found_s1 = searchBoundaryIntersection(self, [t_left, t_right]', [0, NaN]', [NaN, NaN]',...
                                            inlink(1), 'downstream', searchDepth, d_t);
                       
                        t_found_s2 = searchBoundaryIntersection(self, [t_left, t_right]', [0, NaN]', [NaN, NaN]',...
                                            inlink(2), 'downstream', searchDepth, d_t);
                                        
                        t_found_r = searchBoundaryIntersection(self, [t_left, t_right]', [0, NaN]', [NaN, NaN]',...
                                            outlink, 'upstream', searchDepth, d_t);
                        
                        t_found = [t_found_s1; t_found_s2; t_found_r];
                    else
                        t_found = [];
                    end
                    
                    % Aggregate original step intervals if the flows in two
                    % consecutive steps are not changing
                    linkStr = sprintf('link_%d', inlink(1));
                    q_s1 = self.x( self.dv_index.(linkStr).downstream(1,1):...
                                  self.dv_index.(linkStr).downstream(2,1) );
                    linkStr = sprintf('link_%d', inlink(2));
                    q_s2 = self.x( self.dv_index.(linkStr).downstream(1,1):...
                                  self.dv_index.(linkStr).downstream(2,1) );
                    
                    T_tmp_cum_s1 = T_tmp_cum;
                    tmp_group = self.groupSameElement(q_s1, inf);
                    T_tmp_cum_s1 = [ T_tmp_cum_s1(tmp_group(:,1));...
                                     end_time];
                    T_tmp_cum_s2 = T_tmp_cum;
                    tmp_group = self.groupSameElement(q_s2, inf);
                    T_tmp_cum_s2 = [ T_tmp_cum_s2(tmp_group(:,1));...
                                     end_time];
                    
                    T_tmp_cum = unique([T_tmp_cum_s1; T_tmp_cum_s2]);
                                                  
                elseif strcmp( self.net.network_junc.(juncStr).type_junc, 'diverge')
                     
                    if ~isempty(steps.(juncStr))
                        
                        % corner points found from sending link 1 and
                        % receiving link 2, 3
                        t_found_s = searchBoundaryIntersection(self, [t_left, t_right]', [0, NaN]', [NaN, NaN]',...
                                            inlink, 'downstream', searchDepth, d_t);
                        t_found_r1 = searchBoundaryIntersection(self, [t_left, t_right]', [0, NaN]', [NaN, NaN]',...
                                            outlink(1), 'upstream', searchDepth, d_t);
                        t_found_r2 = searchBoundaryIntersection(self, [t_left, t_right]', [0, NaN]', [NaN, NaN]',...
                                            outlink(1), 'upstream', searchDepth, d_t);
                        
                        t_found = [t_found_s; t_found_r1; t_found_r2];
                    
                    else 
                        t_found = [];
                    end
                  
                    % Aggregate original step intervals if the flows in two
                    % consecutive steps are not changing
                    linkStr = sprintf('link_%d', outlink(1));
                    q_r1 = self.x( self.dv_index.(linkStr).upstream(1,1):...
                                  self.dv_index.(linkStr).upstream(2,1) );
                    linkStr = sprintf('link_%d', outlink(2));
                    q_r2 = self.x( self.dv_index.(linkStr).upstream(1,1):...
                                  self.dv_index.(linkStr).upstream(2,1) );
                    
                    T_tmp_cum_r1 = T_tmp_cum;
                    tmp_group = self.groupSameElement(q_r1, inf);
                    T_tmp_cum_r1 = [ T_tmp_cum_r1(tmp_group(:,1));...
                                     end_time];
                    T_tmp_cum_r2 = T_tmp_cum;
                    tmp_group = self.groupSameElement(q_r2, inf);
                    T_tmp_cum_r2 = [ T_tmp_cum_r2(tmp_group(:,1));...
                                     end_time];
                    
                    T_tmp_cum = unique([T_tmp_cum_r1; T_tmp_cum_r2]);
                
                else
                    t_found = [];
                end
                
                % add new found intersection point in the discretization grid
                T_tmp_cum = self.unique_tol([T_tmp_cum; t_found(t_found<end_time)], 1.0e-8);
                
                % return the duration
                T_new_grid.(juncStr).T = T_tmp_cum(2:end) - T_tmp_cum(1:end-1);
                T_new_grid.(juncStr).T_cum = T_tmp_cum;
                
                % Also update the the grid at the connecting links
                for link = self.net.network_junc.(juncStr).inlabel'
                        
                    linkStr = sprintf('link_%d',link);
                    T_new_grid.(linkStr).T_ds = T_tmp_cum(2:end) - T_tmp_cum(1:end-1);
                    T_new_grid.(linkStr).T_ds_cum = T_tmp_cum;
                    T_new_grid.(linkStr).BC_ds = (T_tmp_cum(2:end) - T_tmp_cum(1:end-1))*NaN;
                    
                end
                
                for link = self.net.network_junc.(juncStr).outlabel'
                    
                    linkStr = sprintf('link_%d',link);
                    T_new_grid.(linkStr).T_us = T_tmp_cum(2:end) - T_tmp_cum(1:end-1);
                    T_new_grid.(linkStr).T_us_cum = T_tmp_cum;
                    T_new_grid.(linkStr).BC_us = (T_tmp_cum(2:end) - T_tmp_cum(1:end-1))*NaN;
                    
                end
                
                
            end
            
        end
        
        
        %===============================================================
        function indicedGroupValue = groupSameElement(~, k_ori,minlength)
            % This function is used to group the same consecutive elements in an array into blocks
            % input:
            %       k_ori: a row vector that we would like to group same elements into
            %              a group.
            %       minlength: the minimal length in each group, use inf to disable
            %              this feature.
            % output:
            %       indicedGroupValue: [starting index, end index, block value]
            % example:
            % k=[1 1 1 1 1 1 1 2 2 2 NaN NaN NaN NaN 1 1 1 1 1 2]
            % indicedGroupValue = groupSameElement(k, 5)
            % OUTPUT:
            % indicedGroupValue =
            % [1 5 1;
            %  6 7 1;
            %  8 10 2;
            %  11 14 NaN;
            %  15 19 1
            %  20 20 2];
            
            k_vec = k_ori;
            
            if (~isrow(k_vec) && ~iscolumn(k_vec) )
                error('Check k dimension\n');
            elseif (iscolumn(k_vec))
                k_vec = k_vec';
            end
            
            if (~all(isnan(k_vec)==0))
                randvalue = randn*10;
                while ~all(k_vec~=randvalue)
                    randvalue = randn*10;
                end
                k_vec(isnan(k_vec)) = randvalue;
                randvalue = floor(randvalue*100000)/100000;
            end
            
            %eliminate possible numerical error
            k_tmp = round(k_vec*100000)/100000;
            originalSize = size(k_tmp,2);
            
            i = 1;
            indicedGroupValue = [0 0 0];    %initial offset
            while (1)
                
                i=i+1;
                [values, ivalues, ik_tmp] = unique(k_tmp,'stable');
                
                % set the starting indice
                indicedGroupValue(i,1) = indicedGroupValue(i-1,2)+1;
                % set the value
                indicedGroupValue(i,3) = values(1,1);
                
                if (size(values,2) > 1) %not last
                    if (ivalues(2,1) <= minlength)
                        % in case that we prefer a minlength: e.g. each group
                        % aggregate at most 5 elements.
                        indicedGroupValue(i,2) = ivalues(2,1)-1+ indicedGroupValue(i-1,2);
                        k_tmp(:,1:ivalues(2,1)-1) = [];
                    elseif (ivalues(2,1) > minlength)
                        indicedGroupValue(i,2) = minlength + indicedGroupValue(i-1,2);
                        % remove those processed elements
                        k_tmp(:,1:minlength) = [];
                    end
                    
                    % is the last element, and smaller than minlength
                elseif(originalSize-indicedGroupValue(i,1) <= minlength)
                    indicedGroupValue(i,2) = originalSize;
                    break;
                    
                else
                    indicedGroupValue(i,2) = minlength + indicedGroupValue(i-1,2);
                    % remove those processed elements
                    k_tmp(:,1:minlength) = [];
                end
                
            end
            
            % remove the initial offset (first row)
            indicedGroupValue(1,:) = [];
            
            if (~all(isnan(k_ori)==0))
                indicedGroupValue(indicedGroupValue(:,3)==randvalue,3) = NaN;
            end
            
    
        end
        
        
        %===============================================================
        function abs_e = getAbsError(~, true_sol, approx_sol)
            % This function computes the absolute difference between two solutions.
            % input:
            %       true_sol: struct, .T_cum, .q; the true solution
            %       approx_sol: struct, .T_cum, .q: the true solution
            % output:
            %       abs_e: the absolute error between two solutions
            
            % get combined time grid
            T_grid =[0; unique([true_sol.T_cum; approx_sol.T_cum])];
            
            % compute the abs difference
            abs_e = 0;
            for pt = 2:length(T_grid)
                
                true_flow = true_sol.q( find( true_sol.T_cum>= T_grid(pt),1 ) );
                approx_flow = approx_sol.q( find( approx_sol.T_cum>= T_grid(pt),1 ) );
                abs_e = abs_e + abs( true_flow - approx_flow )*...
                                        (T_grid(pt) - T_grid(pt-1) );
                
            end
            
            
        end
        
        
    end
    
    methods (Access = protected)

        %===============================================================
        function t_found = searchBoundaryIntersection(self, t_interval, M_interval, slope, ...
                                              link, bound, searchDepth, dt_tol)
            % Find the intersection points of the sending and receiving function at the boundary
            %   - this is an older version. The more general function is searchShocksOnLine()                              
            %   - Do NOT call externally.
            % input: 
            %        t_interval, float vector, [t_start; t_end]
            %        M_interval, float vector, [M_start; M_end], NaN if not
            %           computed
            %        link, int, the link id which defines f
            %        bound, 'upstream', 'downstream', which defines f
            %        searchDepth, int, keeps track of the search depth
            %        dr_tol, float, the minimal resolution we will investigate
            % output: found time points 1 x n column
        
            t_found = [];                  
            % if search depth is greater than 3, 
            if searchDepth > self.max_search_depth
                return
            end

            % or the time interval is shorter than dt_res, return middle point  
            if t_interval(2)-t_interval(1) <= dt_tol
                % simply return the middle point as the approximated intersection point
                t_found = sum(t_interval)/2;
                return
            end
            
            % First, compute the two boundary values and the center point
            t_L = t_interval(1);
            t_R = t_interval(2);
            t_C = (t_L+t_R)/2;
            % If value not computed yet
            if isnan(M_interval(1))
                if strcmp(bound, 'downstream')
                    M_L = self.sendingFuncAtLink(t_L, link);
                elseif strcmp(bound, 'upstream')
                    M_L = self.receivingFuncAtLink(t_L, link);
                end
            else
                M_L = M_interval(1);
            end
            if isnan(M_interval(2))
                if strcmp(bound, 'downstream')
                    M_R = self.sendingFuncAtLink(t_R, link);
                elseif strcmp(bound, 'upstream')
                    M_R = self.receivingFuncAtLink(t_R, link);
                end
            else
                M_R = M_interval(2);
            end
            
            if strcmp(bound, 'downstream')
                M_C = self.sendingFuncAtLink(t_C, link);
            elseif strcmp(bound, 'upstream')
                M_C = self.receivingFuncAtLink(t_C, link);
            end
            
            % Second, check if on a straight line, return []
            if self.onStraightLine([t_L; t_C; t_R], [M_L; M_C; M_R])
                % double check if on a straight line. Some times, the
                % segments are symmetric and those 3 points are on the same
                % line, but the sending function is not a line
                t_Q_1 = t_L + (t_R-t_L)/4;      % first quarter
                t_Q_3 = t_L + 3*(t_R-t_L)/4;    % third quarter
                
                if strcmp(bound, 'downstream')
                    M_Q_1 = self.sendingFuncAtLink(t_Q_1, link);
                    M_Q_3 = self.sendingFuncAtLink(t_Q_3, link);
                elseif strcmp(bound, 'upstream')
                    M_Q_1 = self.receivingFuncAtLink(t_Q_1, link);
                    M_Q_3 = self.receivingFuncAtLink(t_Q_3, link);
                end

                if self.onStraightLine([t_L; t_Q_1; t_C], [M_L; M_Q_1; M_C]) &&...
                   self.onStraightLine([t_C; t_Q_3; t_R], [M_C; M_Q_3; M_R])
                    % sampled 5 points on the same line, then very likely
                    % to be a straight line
                    return
                end
            end
            
            % Third, find the left and right slope and compute the
            % intersection 
            if isnan(slope(1))
                s_L = self.findBoundaryFuncSlope([t_L; t_C], [M_L; M_C], 'left',...
                                 link, bound, dt_tol);
            else
                s_L = slope(1);
            end
            if isnan(slope(2))
                s_R = self.findBoundaryFuncSlope([t_C; t_R], [M_C; M_R], 'right',...
                            	 link, bound, dt_tol);
            else
                s_R = slope(2);
            end
            % intersection point 
            if abs(s_L - s_R) >= 1.0e-3
                % make sure two slopes are not parallel, otherwise directly
                % split them.
                t_insct = (M_R - M_L + s_L*t_L - s_R*t_R)/(s_L - s_R);
                M_insct = s_L*(t_insct - t_L) + M_L;
                
                if t_insct <= t_R && t_insct >= t_L
                    % if the intersection is in the time step, then
                    % validate, otherwise, move on to split
                    if strcmp(bound, 'downstream')
                        M_validate = self.sendingFuncAtLink(t_insct, link);
                    elseif strcmp(bound, 'upstream')
                        M_validate = self.receivingFuncAtLink(t_insct, link);
                    end
                    
                    if abs(M_insct - M_validate) <= 1e-6
                        % found the intersection
                        t_found = [t_found; t_insct];
                        return
                    end
                end
            end
            
            % Forth, recursion, split the interval if there are more than 
            % one corners based on: 
            % 1. slopes are parallel; 2. intersection outside of
            % the interval; 3. cornerpoint not valid in sample
            t_found_left = self.searchBoundaryIntersection([t_L; t_C], [M_L; M_C],...
                [NaN; NaN], link, bound, searchDepth+1, dt_tol);
            t_found_right = self.searchBoundaryIntersection([t_C; t_R], [M_C; M_R],...
                [NaN; NaN], link, bound, searchDepth+1, dt_tol);
            
            t_found = [t_found; t_found_left; t_found_right];
                        
        end
        
        
        %===============================================================
        function slope = findBoundaryFuncSlope(self, t_interval, M_interval, slope_side, ...
                                   link, bound, dt_tol)
            % Find the left or right slope at the boundary point of a piecewise linear line
            %   - Do NOT call externally
            % input: 
            %        t_interval, float vector, [t_start; t_end]
            %        M_interval, float vector, [M_start; M_end], NaN if not
            %           computed
            %        slope_side, 'left','right', which side slope
            %        link, int, the link id which defines f
            %        bound, 'upstream', 'downstream', which defines f
            %        dr_tol, float, the minimal resolution we will investigate
            % output: found time points 1 x n column
                               
            % if the time interval is too small, than stop computing slope
            % since this may give large numerical error
            if t_interval(2) - t_interval(1) <= dt_tol
                slope = NaN;
                return
            end
            
            % First, sample boundary points values and the center point
            t_L = t_interval(1);
            t_R = t_interval(2);
            t_C = (t_L+t_R)/2;
            % If value not computed yet
            if isnan(M_interval(1))
                if strcmp(bound, 'downstream')
                    M_L = self.sendingFuncAtLink(t_L, link);
                elseif strcmp(bound, 'upstream')
                    M_L = self.receivingFuncAtLink(t_L, link);
                end
            else
                M_L = M_interval(1);
            end
            if isnan(M_interval(2))
                if strcmp(bound, 'downstream')
                    M_R = self.sendingFuncAtLink(t_R, link);
                elseif strcmp(bound, 'upstream')
                    M_R = self.receivingFuncAtLink(t_R, link);
                end
            else
                M_R = M_interval(2);
            end
            
            if strcmp(bound, 'downstream')
                M_C = self.sendingFuncAtLink(t_C, link);
            elseif strcmp(bound, 'upstream')
                M_C = self.receivingFuncAtLink(t_C, link);
            end
            
            % Second, check if three points are on the same line
            if self.onStraightLine([t_L; t_C; t_R],[M_L; M_C; M_R])
                % if true, return the slope
                slope = (M_R-M_L)/(t_R-t_L);
            else
                if strcmp(slope_side, 'left')
                    % keep the left half to find the slope
                    slope = self.findBoundaryFuncSlope([t_L; t_C], [M_L; M_C], slope_side,...
                                         link, bound, dt_tol);
                elseif strcmp(slope_side, 'right')
                    % keep the right half to find the slope
                    slope = self.findBoundaryFuncSlope([t_C; t_R], [M_C; M_R], slope_side,...
                                         link, bound, dt_tol);
                end
            end                                 
        end
        
        
        
        %===============================================================
        function TF = onStraightLine(~, x, y)
            % This function tells whether three points are on the same straight line
            % input:
            %       x: 3x1 vector, float
            %       y: 3x1 vector, float
            % output:
            %       TF: true if three points on one straight line with small
            %           tolerance; otherwise false
            
            % if a horizontal line
            if abs(x(3)-x(2)) <= 1.0e-10 && abs(x(2)-x(1)) <= 1.0e-6
                TF = true;
                return
            end
            
            % if a vertical line
            if abs(y(3)-y(2)) <= 1.0e-10 && abs(y(2)-y(1)) <= 1.0e-6
                TF = true;
                return
            end
            
            % otherwise, compute slope
            slope = (y(3)-y(1))/(x(3)-x(1));
            if abs(slope*(x(2)-x(1)) + y(1) - y(2)) <= 1.0e-6
                TF = true;
            else
                TF = false;
            end
            
        end
        
        
        
        %===============================================================
        function d_M = entropicSolutionAtJunc(self, t_sample, junc)
            % This function computes the admissible solution between the self.t_ref and t_sample
            %   - it samples the sending and receiving function at all links,
            %   - then it computes the admissible solution
            %   - do NOT use externally
            % input: 
            %        t_sample: float, the time point to be sampled
            %        junc: the label of the junction to be sampled
            % output: 
            %       d_M: for connection, float, the entropic car ID at t_sample
            %            for merge/diverge, 2x1 float, the entropic car ID at
            %            two incoming/outgoing links, with respect to self.t_ref
            
            % t_ref is the reference point where M = 0
            if t_sample == self.t_ref
                d_M = 0;
                return
            end
            
            juncStr = sprintf('junc_%d',junc);
            
            if strcmp(self.net.network_junc.(juncStr).type_junc,'connection')
                
                inLink = self.net.network_junc.(juncStr).inlabel;
                outLink = self.net.network_junc.(juncStr).outlabel;
                
                % simply take the minimum of the sending and receiving of
                % two links
                
                d_M = min(self.sendingFuncAtLink(t_sample, inLink),...
                          self.receivingFuncAtLink(t_sample, outLink));
                      
            % if it is a onrampjunc, then the entropic solution should be
            % q_upFreeway = min(q_upSending, q_downReceiving - q_on)
            elseif strcmp(self.net.network_junc.(juncStr).type_junc, 'onrampjunc')
                
                % find the upstream freeway, onramp, and the downstream
                % freeway
                inlinks = self.net.network_junc.(juncStr).inlabel;
                outlink = self.net.network_junc.(juncStr).outlabel;
                linkStr1 = sprintf('link_%d',inlinks(1));
                linkStr2 = sprintf('link_%d',inlinks(2));
                
                if strcmp(self.net.network_hwy.(linkStr1).para_linktype, 'freeway')
                    onrampStr = linkStr2;
                    up_link = inlinks(1);
                else
                    onrampStr = linkStr1;
                    up_link = inlinks(2);
                end
                
                % extract the controlled onramp flow at the investigated
                % step end at time t_sample
                q_on = self.x(self.dv_index.(onrampStr).downstream(1,1):...
                                   self.dv_index.(onrampStr).downstream(2,1));
                T_cum_grid = self.net.network_junc.(juncStr).T_cum;
                T_grid = self.net.network_junc.(juncStr).T;
                step = find(T_cum_grid == t_sample) - 1;
                               
                d_M = min( self.sendingFuncAtLink(t_smaple, up_link), ...
                           self.receivingFuncAtLink(t_sample, outlink) -...
                                q_on(step)*T_grid(step));
                            
                            
            % if it is a offrampjunc, then the entropic solution should be
            % q_downFreeway = min(q_downReceiving, q_upSending - q_off)
            elseif strcmp(self.net.network_junc.(juncStr).type_junc, 'offrampjunc')
                
                % find the upstream freeway, onramp, and the downstream
                % freeway
                inlink = self.net.network_junc.(juncStr).inlabel;
                outlinks = self.net.network_junc.(juncStr).outlabel;
                linkStr1 = sprintf('link_%d',outlinks(1));
                linkStr2 = sprintf('link_%d',outlinks(2));
                
                if strcmp(self.net.network_hwy.(linkStr1).para_linktype, 'freeway')
                    offrampStr = linkStr2;
                    down_link = outlinks(1);
                else
                    offrampStr = linkStr1;
                    down_link = outlinks(2);
                end
                
                % extract the controlled offramp flow at the investigated
                % step end at time t_sample
                q_off = self.x(self.dv_index.(offrampStr).upstream(1,1):...
                                   self.dv_index.(offrampStr).upstream(2,1));
                T_cum_grid = self.net.network_junc.(juncStr).T_cum;
                T_grid = self.net.network_junc.(juncStr).T;
                step = find(T_cum_grid == t_sample) - 1;
                               
                d_M = min( self.receivingFuncAtLink(t_smaple, down_link), ...
                           self.sendingFuncAtLink(t_sample, inlink) -...
                                q_off(step)*T_grid(step));    
            
            elseif strcmp(self.net.network_junc.(juncStr).type_junc,'merge')
                
                inLinks = self.net.network_junc.(juncStr).inlabel;
                outLink = self.net.network_junc.(juncStr).outlabel;
                
                % slope of the priority with y: inLinks(2); and x:
                % inlinks(1)
                sPriority = self.net.network_junc.(juncStr).ratio(2)/...
                    self.net.network_junc.(juncStr).ratio(1);
                
                % here we extract the sending and receiving function
                M_R = self.receivingFuncAtLink(t_sample, outLink);
                M_S1 = self.sendingFuncAtLink(t_sample, inLinks(1));
                M_S2 = self.sendingFuncAtLink(t_sample, inLinks(2));
                
                % case 1: M_S1 + M_S2 <= M_R, then through flow is
                % M_s1 + M_s2
                if M_S1 + M_S2 <= M_R
                    d_M = [M_S1; M_S2];
                
                % case 2: intersection point in side (0-M_S1, 0-M_S2) box
                elseif ( M_R/(1+sPriority) <= M_S1) &&...
                        (sPriority*M_R/(1+sPriority) <= M_S2)
                    
                    d_M = [M_R/(1+sPriority); sPriority*M_R/(1+sPriority)];                    
                
                % case 3: constrained by M_S1 flow
                elseif M_R/(1+sPriority) > M_S1
                    d_M = [M_S1, M_R-M_S1];
                
                % case 4: constrained by M_S2 flow
                elseif sPriority*M_R/(1+sPriority) > M_S2
                    d_M = [M_R-M_S2 , M_S2];
                end

            elseif strcmp(self.net.network_junc.(juncStr).type_junc,'diverge')
                
                inLink = self.net.network_junc.(juncStr).inlabel;
                outLinks = self.net.network_junc.(juncStr).outlabel;
                
                % slope of the distribution with y: outLinks(2); and x:
                % outlinks(1)
                sDistribution = self.net.network_junc.(juncStr).ratio(2)/...
                    self.net.network_junc.(juncStr).ratio(1);
                
                % here we extract the sending and receiving function
                M_R1 = self.receivingFuncAtLink(t_sample, outLinks(1));
                M_R2 = self.receivingFuncAtLink(t_sample, outLinks(2));
                M_S = self.sendingFuncAtLink(t_samplw, inLink);
                
                % case 1: M_R1 + M_R2 <= M_S, then through flow is
                % M_R1+M_R2
                if M_R1 + M_R2 <= M_S
                    d_M = [M_R1, M_R2];
                    
                % case 2: intersection point in side (0-M_R1, 0-M_R2) box
                elseif ( M_S/(1+sDistribution) <= M_R1) &&...
                        (sDistribution*M_S/(1+sDistribution) <= M_R2)
                    d_M = [M_S/(1+sDistribution), sDistribution*M_S/(1+sDistribution)];
                    
                % case 3: constrained by M_R1 flow
                elseif M_S/(1+sDistribution) > M_R1
                    d_M = [M_R1, M_S-M_R1];
                    
                % case 4: constrained by M_R2 flow
                elseif sDistribution*M_S/(1+sDistribution) > M_R2
                    d_M = [M_S-M_R2 , M_R2];
                end

            end
            
        end
        
        
        
        %===============================================================
        function d_M = sendingFuncAtLink(self, t_sample, link)
            % This function cmoputes the amount of vehicle that can be sent on link between time self.t_ref and t_sample
            %   - The closest boundary condition (solution from x) was removed 
            %       since it may not be entropic.
            %   - Do NOT call externally.
            % input: 
            %        t_sample: float, the boundary time to be sampled
            %        link: int, the link ID to be sampled
            % output: d_M: the amount of vehicles that can be sent
            
            % since t_ref is the reference point, so simply 0
            if t_sample == self.t_ref
                d_M = 0;
                return
            end
            
            % save the absolute vehicle id on this link in M = [ 0; 0 ]
            M = ones(2,1)*NaN;
            
            % compute the M of t on this link
            %===============================================================
            % extract the fundamental diagram
            linkStr = sprintf('link_%d',link);
            v_f = self.net.network_hwy.(linkStr).para_vf;
            w = self.net.network_hwy.(linkStr).para_w;    % < 0
            k_c = self.net.network_hwy.(linkStr).para_kc;
            k_m = self.net.network_hwy.(linkStr).para_km;            
            
            % Get the initial number of vehicles 
            IC(:,1) = self.net.network_hwy.(linkStr).X_grid_cum(1:end-1);
            IC(:,2) = self.net.network_hwy.(linkStr).X_grid_cum(2:end);
            IC(:,3) = self.net.network_hwy.(linkStr).IC;
            IC_num_veh = - ( IC(:,2) - IC(:,1)).*IC(:,3);
            IC_cum_num_veh = [0; cumsum(IC_num_veh) ];
            IC_total_num_veh = IC_cum_num_veh(end); %  should be negative
            
            % remove the last entry, which overlaps with the downstream
            % conditions
            IC_cum_num_veh(end) = [];
   
            if ~iscolumn(IC_cum_num_veh)
                IC_cum_num_veh = IC_cum_num_veh';
            end
            
            % select the boudnary conditions we need to use to compute the
            % vehicle id at the sampling point.
            % We only need to use the boundary conditions both upstream and
            % downstream up to the sampling point time.
            % the number of boundary steps we need to take into account
            num_us_pre_steps = sum(self.net.network_hwy.(linkStr).T_us_cum < t_sample);
            num_ds_pre_steps = sum(self.net.network_hwy.(linkStr).T_ds_cum < t_sample);
            
            
            % First extract the boundary conditions that to be considered
            % for computing the vehilce label at the sampling point
            BC_us(:,1) = self.net.network_hwy.(linkStr).T_us_cum(1:num_us_pre_steps);
            BC_us(:,2) = self.net.network_hwy.(linkStr).T_us_cum(2:num_us_pre_steps+1);
            BC_us(num_us_pre_steps, 2) = t_sample;
            BC_ds(:,1) = self.net.network_hwy.(linkStr).T_ds_cum(1:num_ds_pre_steps);
            BC_ds(:,2) = self.net.network_hwy.(linkStr).T_ds_cum(2:num_ds_pre_steps+1);
            BC_ds(num_ds_pre_steps, 2) = t_sample;
            
            % Use the computed solution values as the boundary condition
            BC_us(:,3) = self.x( self.dv_index.(linkStr).upstream(1,1):...
                                self.dv_index.(linkStr).upstream(1,1) + num_us_pre_steps-1);
            BC_ds(:,3) = self.x( self.dv_index.(linkStr).downstream(1,1):...
                                self.dv_index.(linkStr).downstream(1,1) + num_ds_pre_steps-1);
            
            % Second compute the absolute vehicle ID
            BC_us_num_veh = (BC_us(:,2) - BC_us(:,1)).*BC_us(:,3);
            BC_cum_us_M = [0; cumsum(BC_us_num_veh)];
            BC_cum_us_M(end) = [];
            
            BC_ds_num_veh = (BC_ds(:,2) - BC_ds(:,1)).*BC_ds(:,3);
            BC_cum_ds_M = [0; cumsum(BC_ds_num_veh) ] + IC_total_num_veh;
            BC_cum_ds_M(end) = [];
            
            % remove the nonentropic boundary condition
            if sum( self.net.network_hwy.(linkStr).T_ds_cum > self.t_ref &...
                    self.net.network_hwy.(linkStr).T_ds_cum < t_sample ) ~= 0
                % remove the last two pieces of downstream conditions
                BC_ds( max(1,num_ds_pre_steps-1):num_ds_pre_steps , 3) = NaN;
            else
                % remove the last piece of downstream condition
                BC_ds( num_ds_pre_steps , 3) = NaN;
            end
            
            % points coordinates
            len_link = self.net.network_hwy.(linkStr).para_postkm*1000;
            position = len_link;
            
            t_array = [self.t_ref; t_sample];
            for i = 1:length(t_array)
                %==========================================================
                % compute the solution associated to the initial condition
                lowerThanKc = (IC(:,3) <= k_c);
                greaterThanKc = (IC(:,3) > k_c);
                
                %INITIAL==== <k_c ===========================================
                % first check those conditions with <= k_c
                inCharacterDomain = IC(:,1)+v_f*t_array(i) <= position &...
                                    IC(:,2)+v_f*t_array(i) >= position;
                
                % only applies to those <k_c
                activeCharacterDomain = lowerThanKc & inCharacterDomain;
                
                inFanDomain = IC(:,1)+v_f*t_array(i) > position &...
                              IC(:,1)+w*t_array(i) <= position;
                
                activeFanDomain = lowerThanKc & inFanDomain;
                
                % solution in characteristic domain
                tmp_M = min( IC_cum_num_veh(activeCharacterDomain,1) -...
                    IC(activeCharacterDomain,3).*(position - v_f*t_array(i) - ...
                    IC(activeCharacterDomain,1)) );
                M(i) = self.minNonEmpty(M(i), tmp_M);
                % solution in fan domain
                tmp_M = min( IC_cum_num_veh(activeFanDomain,1) -...
                    k_c.*(position - v_f*t_array(i) - ...
                    IC(activeFanDomain,1)  )  );
                M(i) = self.minNonEmpty(M(i), tmp_M);
                
                %INITIAL==== >k_c ===========================================
                % now check the initial condition with > k_c
                inCharacterDomain = IC(:,1)+w*t_array(i) <= position &...
                                    IC(:,2)+w*t_array(i) >= position; 
                activeCharacterDomain = greaterThanKc & inCharacterDomain;
                
                inFanDomain = IC(:,2)+v_f*t_array(i) > position &...
                              IC(:,2)+w*t_array(i) <= position; 
                activeFanDomain = greaterThanKc & inFanDomain;
                
                % solution in characteristic domain
                tmp_M = min( IC_cum_num_veh(activeCharacterDomain,1) -...
                    IC(activeCharacterDomain,3).*(position - w*t_array(i) - ...
                    IC(activeCharacterDomain,1)) -k_m*t_array(i)*w);
                M(i) = self.minNonEmpty(M(i), tmp_M);
                % solution in fan domain
                tmp_M = min( IC_cum_num_veh(activeFanDomain,1)+ IC_num_veh(activeFanDomain,1)...
                    - k_c.*(position - w*t_array(i) - IC(activeFanDomain,2)) ...
                    - k_m*t_array(i)*w);
                M(i) = self.minNonEmpty(M(i), tmp_M);
                
                %UPSTREAM==================================================
                % now check the upstream conditions, only apply those
                % boundary conditions that are not NaN
                notNaN = ~isnan(BC_us(:,3));
                
                inCharacterDomain = BC_us(:,1)+ position/v_f <= t_array(i) &...
                                    BC_us(:,2)+ position/v_f >= t_array(i); 
                activeCharacterDomain = notNaN & inCharacterDomain;
                
                inFanDomain = BC_us(:,2)+ position/v_f < t_array(i);
                activeFanDomain = notNaN & inFanDomain;
                
                % solution in characteristic domain
                tmp_M = min( BC_cum_us_M(activeCharacterDomain,1) +...
                    BC_us(activeCharacterDomain,3).*(t_array(i) - position/v_f - ...
                    BC_us(activeCharacterDomain,1)) );
                M(i) = self.minNonEmpty(M(i), tmp_M);
                % solution in fan domain
                tmp_M = min( BC_cum_us_M(activeFanDomain,1) +...
                    BC_us_num_veh(activeFanDomain,1) + ...
                    k_c*v_f.*(t_array(i) - position/v_f - ...
                    BC_us(activeFanDomain,2)) );
                M(i) = self.minNonEmpty(M(i), tmp_M);
                
                
                %DOWNSTREAM==================================================
                % now check the upstream conditions
                notNaN = ~isnan(BC_ds(:,3));
                
                inCharacterDomain = BC_ds(:,1) + (position - len_link)/w <= t_array(i) &...
                                    BC_ds(:,2)+ (position - len_link)/w >= t_array(i);
                activeCharacterDomain = notNaN & inCharacterDomain;
                
                inFanDomain = BC_ds(:,2)+ (position - len_link)/w < t_array(i);
                activeFanDomain = notNaN & inFanDomain;
                
                % solution in characteristic domain
                tmp_M = min( BC_cum_ds_M(activeCharacterDomain,1) +...
                    BC_ds(activeCharacterDomain,3).*(t_array(i) -...
                    (position - len_link)/w - ...
                    BC_ds(activeCharacterDomain,1))...
                    - k_m*(position - len_link));
                M(i) = self.minNonEmpty(M(i), tmp_M);
                % solution in fan domain
                tmp_M = min( BC_cum_ds_M(activeFanDomain,1) +...
                    BC_ds_num_veh(activeFanDomain,1) + ...
                    k_c*v_f.*(t_array(i) - (position-len_link)/v_f - ...
                    BC_ds(activeFanDomain,2)) );
                
                M(i) = self.minNonEmpty(M(i), tmp_M);
                
                
            end % end for each t (here we just have (t_ref and t) )
            
            if all(~isnan(M)) && all(~isempty(M))
                
                % use car ID at t_ref as a reference
                d_M = M(2) - M(1);
                
            else
                error('failed to sample points');
            end
            
            
        end 
        
        
        
        %===============================================================
        function d_M = receivingFuncAtLink(self, t_sample, link)
            % This function cmoputes the amount of vehicle that can be received on link between time self.t_ref and t_sample
            %   - The closest boundary condition (solution from x) was removed 
            %       since it may not be entropic.
            %   - Do NOT call externally.
            % input: 
            %        t_sample: float, the boundary time to be sampled
            %        link: int, the link ID to be sampled
            % output: d_M, the amount of vehicles that can be received
            
            % since t_ref is the reference point, so simply 0
            if t_sample == self.t_ref
                d_M = 0;
                return
            end
            
            % save the absolute vehicle id on this link in M = [ 0; 0 ]
            M = ones(2,1)*NaN;
            
            % compute the M of t on this link
            %===============================================================
            % extract the fundamental diagram
            linkStr = sprintf('link_%d',link);
            v_f = self.net.network_hwy.(linkStr).para_vf;
            w = self.net.network_hwy.(linkStr).para_w;    % < 0
            k_c = self.net.network_hwy.(linkStr).para_kc;
            k_m = self.net.network_hwy.(linkStr).para_km;            
            
            % Get the initial number of vehicles 
            IC(:,1) = self.net.network_hwy.(linkStr).X_grid_cum(1:end-1);
            IC(:,2) = self.net.network_hwy.(linkStr).X_grid_cum(2:end);
            IC(:,3) = self.net.network_hwy.(linkStr).IC;
            IC_num_veh = - ( IC(:,2) - IC(:,1)).*IC(:,3);
            IC_cum_num_veh = [0; cumsum(IC_num_veh) ];
            IC_total_num_veh = IC_cum_num_veh(end); %  should be negative
            
            % remove the last entry, which overlaps with the downstream
            % conditions
            IC_cum_num_veh(end) = [];
   
            if ~iscolumn(IC_cum_num_veh)
                IC_cum_num_veh = IC_cum_num_veh';
            end
            
            % select the boudnary conditions we need to use to compute the
            % vehicle id at the sampling point.
            % We only need to use the boundary conditions both upstream and
            % downstream up to the sampling point time.
            % the number of boundary steps we need to take into account
            num_us_pre_steps = sum(self.net.network_hwy.(linkStr).T_us_cum < t_sample);
            num_ds_pre_steps = sum(self.net.network_hwy.(linkStr).T_ds_cum < t_sample);
            
            
            % First extract the boundary conditions that to be considered
            % for computing the vehilce label at the sampling point
            BC_us(:,1) = self.net.network_hwy.(linkStr).T_us_cum(1:num_us_pre_steps);
            BC_us(:,2) = self.net.network_hwy.(linkStr).T_us_cum(2:num_us_pre_steps+1);
            BC_us(num_us_pre_steps, 2) = t_sample;
            BC_ds(:,1) = self.net.network_hwy.(linkStr).T_ds_cum(1:num_ds_pre_steps);
            BC_ds(:,2) = self.net.network_hwy.(linkStr).T_ds_cum(2:num_ds_pre_steps+1);
            BC_ds(num_ds_pre_steps, 2) = t_sample;
            
            % Use the computed solution values as the boundary condition
            BC_us(:,3) = self.x( self.dv_index.(linkStr).upstream(1,1):...
                                self.dv_index.(linkStr).upstream(1,1) + num_us_pre_steps-1);
            BC_ds(:,3) = self.x( self.dv_index.(linkStr).downstream(1,1):...
                                self.dv_index.(linkStr).downstream(1,1) + num_ds_pre_steps-1);
            
            % Second compute the absolute vehicle ID
            BC_us_num_veh = (BC_us(:,2) - BC_us(:,1)).*BC_us(:,3);
            BC_cum_us_M = [0; cumsum(BC_us_num_veh)];
            BC_cum_us_M(end) = [];
            
            BC_ds_num_veh = (BC_ds(:,2) - BC_ds(:,1)).*BC_ds(:,3);
            BC_cum_ds_M = [0; cumsum(BC_ds_num_veh) ] + IC_total_num_veh;
            BC_cum_ds_M(end) = [];
            
            % remove nonentropic upstream boundary condition
            if sum( self.net.network_hwy.(linkStr).T_us_cum > self.t_ref &...
                    self.net.network_hwy.(linkStr).T_us_cum < t_sample) ~= 0
                % remove the last two pieces of upstream conditions
                BC_us( max(1, num_us_pre_steps-1) : num_us_pre_steps, 3) = NaN;
            else
                % remove the last piece of upstream condition
                BC_us( num_us_pre_steps, 3) = NaN;
            end
            
            len_link = self.net.network_hwy.(linkStr).para_postkm*1000;
            position = 0;
           
            t_array = [self.t_ref; t_sample];
            for i = 1:length(t_array)
                %==========================================================
                % compute the solution associated to the initial condition
                lowerThanKc = (IC(:,3) <= k_c);
                greaterThanKc = (IC(:,3) > k_c);
                
                %INITIAL==== <k_c ===========================================
                % first check those conditions with <= k_c
                inCharacterDomain = IC(:,1)+v_f*t_array(i) <= position &...
                                    IC(:,2)+v_f*t_array(i) >= position;
                
                % only applies to those <k_c
                activeCharacterDomain = lowerThanKc & inCharacterDomain;
                
                inFanDomain = IC(:,1)+v_f*t_array(i) > position &...
                              IC(:,1)+w*t_array(i) <= position;
                
                activeFanDomain = lowerThanKc & inFanDomain;
                
                % solution in characteristic domain
                tmp_M = min( IC_cum_num_veh(activeCharacterDomain,1) -...
                    IC(activeCharacterDomain,3).*(position - v_f*t_array(i) - ...
                    IC(activeCharacterDomain,1)) );
                M(i) = self.minNonEmpty(M(i), tmp_M);
                % solution in fan domain
                tmp_M = min( IC_cum_num_veh(activeFanDomain,1) -...
                    k_c.*(position - v_f*t_array(i) - ...
                    IC(activeFanDomain,1)  )  );
                M(i) = self.minNonEmpty(M(i), tmp_M);
                
                %INITIAL==== >k_c ===========================================
                % now check the initial condition with > k_c
                inCharacterDomain = IC(:,1)+w*t_array(i) <= position &...
                                    IC(:,2)+w*t_array(i) >= position; 
                activeCharacterDomain = greaterThanKc & inCharacterDomain;
                
                inFanDomain = IC(:,2)+v_f*t_array(i) > position &...
                              IC(:,2)+w*t_array(i) <= position; 
                activeFanDomain = greaterThanKc & inFanDomain;
                
                % solution in characteristic domain
                tmp_M = min( IC_cum_num_veh(activeCharacterDomain,1) -...
                    IC(activeCharacterDomain,3).*(position - w*t_array(i) - ...
                    IC(activeCharacterDomain,1)) -k_m*t_array(i)*w);
                M(i) = self.minNonEmpty(M(i), tmp_M);
                % solution in fan domain
                tmp_M = min( IC_cum_num_veh(activeFanDomain,1)+ IC_num_veh(activeFanDomain,1)...
                    - k_c.*(position - w*t_array(i) - IC(activeFanDomain,2)) ...
                    - k_m*t_array(i)*w);
                M(i) = self.minNonEmpty(M(i), tmp_M);
                
                %UPSTREAM==================================================
                % now check the upstream conditions, only apply those
                % boundary conditions that are not NaN
                notNaN = ~isnan(BC_us(:,3));
                
                inCharacterDomain = BC_us(:,1)+ position/v_f <= t_array(i) &...
                                    BC_us(:,2)+ position/v_f >= t_array(i); 
                activeCharacterDomain = notNaN & inCharacterDomain;
                
                inFanDomain = BC_us(:,2)+ position/v_f < t_array(i);
                activeFanDomain = notNaN & inFanDomain;
                
                % solution in characteristic domain
                tmp_M = min( BC_cum_us_M(activeCharacterDomain,1) +...
                    BC_us(activeCharacterDomain,3).*(t_array(i) - position/v_f - ...
                    BC_us(activeCharacterDomain,1)) );
                M(i) = self.minNonEmpty(M(i), tmp_M);
                % solution in fan domain
                tmp_M = min( BC_cum_us_M(activeFanDomain,1) +...
                    BC_us_num_veh(activeFanDomain,1) + ...
                    k_c*v_f.*(t_array(i) - position/v_f - ...
                    BC_us(activeFanDomain,2)) );
                M(i) = self.minNonEmpty(M(i), tmp_M);
                
                
                %DOWNSTREAM==================================================
                % now check the upstream conditions
                notNaN = ~isnan(BC_ds(:,3));
                
                inCharacterDomain = BC_ds(:,1) + (position - len_link)/w <= t_array(i) &...
                                    BC_ds(:,2)+ (position - len_link)/w >= t_array(i);
                activeCharacterDomain = notNaN & inCharacterDomain;
                
                inFanDomain = BC_ds(:,2)+ (position - len_link)/w < t_array(i);
                activeFanDomain = notNaN & inFanDomain;
                
                % solution in characteristic domain
                tmp_M = min( BC_cum_ds_M(activeCharacterDomain,1) +...
                    BC_ds(activeCharacterDomain,3).*(t_array(i) -...
                    (position - len_link)/w - ...
                    BC_ds(activeCharacterDomain,1))...
                    - k_m*(position - len_link));
                M(i) = self.minNonEmpty(M(i), tmp_M);
                % solution in fan domain
                tmp_M = min( BC_cum_ds_M(activeFanDomain,1) +...
                    BC_ds_num_veh(activeFanDomain,1) + ...
                    k_c*v_f.*(t_array(i) - (position-len_link)/v_f - ...
                    BC_ds(activeFanDomain,2)) );
                
                M(i) = self.minNonEmpty(M(i), tmp_M);
                
                
            end % end for each t (here we just have (t_ref and t) )
            
            if all(~isnan(M)) && all(~isempty(M))
                
                % use car ID at t_ref as a reference
                d_M = M(2) - M(1);
                
            else
                error('failed to sample points');
            end
            
            
        end 
        
        
        
        %===============================================================
        function M = absVehID(self, link, times, positions)
            % This function computes the Moskowitz solution at (t,x) given initial, upstream and downstream boudnary conditions.
            % input:
            %       link: the link id to be sampled
            %       time: column vector, the time of the sampled points
            %       position: column, the position of the point in meters
            % output:
            %       M: the absolute vehicle id at the position with M=0 at left bottom
            
            % it must be points
            if length(times) ~= length(positions)
                error('ERROR: the time and position must have the same length\n')
            end

            % save the absolute vehicle id on this link in M = [ 0; 0 ]
            M = ones(length(times),1)*NaN;            
            
            %===============================================================
            % extract the fundamental diagram
            linkStr = sprintf('link_%d',link);
            v_f = self.net.network_hwy.(linkStr).para_vf;
            w = self.net.network_hwy.(linkStr).para_w;    % < 0
            k_c = self.net.network_hwy.(linkStr).para_kc;
            k_m = self.net.network_hwy.(linkStr).para_km;     
            len_link = self.net.network_hwy.(linkStr).para_postkm*1000;
            
            % Get the initial number of vehicles 
            IC(:,1) = self.net.network_hwy.(linkStr).X_grid_cum(1:end-1);
            IC(:,2) = self.net.network_hwy.(linkStr).X_grid_cum(2:end);
            IC(:,3) = self.x( self.dv_index.(linkStr).initial(1,1):...
                              self.dv_index.(linkStr).initial(2,1));
            IC_num_veh = - ( IC(:,2) - IC(:,1)).*IC(:,3);
            IC_cum_num_veh = [0; cumsum(IC_num_veh) ];
            IC_total_num_veh = IC_cum_num_veh(end); %  should be negative
            
            % remove the last entry, which overlaps with the downstream
            % conditions
            IC_cum_num_veh(end) = [];
   
            if ~iscolumn(IC_cum_num_veh)
                IC_cum_num_veh = IC_cum_num_veh';
            end

            % sample each point
            for i = 1:length(times)

                time = times(i);
                position = positions(i);

                % select the boudnary conditions we need to use to compute the
                % vehicle id at the sampling point.
                % We only need to use the boundary conditions both upstream and
                % downstream up to the sampling point time.
                % the number of boundary steps we need to take into account
                num_us_pre_steps = sum(self.net.network_hwy.(linkStr).T_us_cum < time);
                num_ds_pre_steps = sum(self.net.network_hwy.(linkStr).T_ds_cum < time);

                 % First extract the boundary conditions that to be considered
                % for computing the vehilce label at the sampling point
                BC_us(:,1) = self.net.network_hwy.(linkStr).T_us_cum(1:num_us_pre_steps);
                BC_us(:,2) = self.net.network_hwy.(linkStr).T_us_cum(2:num_us_pre_steps+1);
                BC_us(num_us_pre_steps, 2) = time;
                BC_ds(:,1) = self.net.network_hwy.(linkStr).T_ds_cum(1:num_ds_pre_steps);
                BC_ds(:,2) = self.net.network_hwy.(linkStr).T_ds_cum(2:num_ds_pre_steps+1);
                BC_ds(num_ds_pre_steps, 2) = time;
                
                % Use the computed solution values as the boundary condition
                BC_us(:,3) = self.x( self.dv_index.(linkStr).upstream(1,1):...
                                    self.dv_index.(linkStr).upstream(1,1) + num_us_pre_steps-1);
                BC_ds(:,3) = self.x( self.dv_index.(linkStr).downstream(1,1):...
                                    self.dv_index.(linkStr).downstream(1,1) + num_ds_pre_steps-1);

                % Second compute the absolute vehicle ID at boundaries
                BC_us_num_veh = (BC_us(:,2) - BC_us(:,1)).*BC_us(:,3);
                BC_cum_us_M = [0; cumsum(BC_us_num_veh)];
                BC_cum_us_M(end) = [];
                
                BC_ds_num_veh = (BC_ds(:,2) - BC_ds(:,1)).*BC_ds(:,3);
                BC_cum_ds_M = [0; cumsum(BC_ds_num_veh) ] + IC_total_num_veh;
                BC_cum_ds_M(end) = [];

                % NOTE: Do NOT remove any boundary condition before time.
                % In sendingFuncAtLink(), we remove the closest piece of boundary condition 
                % is because we would like to compute the sending function and the closest 
                % piece of BC is nonentropic. Here we only want to know the vehicle id
                % given the solutoin (which should already be entropic).


                %==========================================================
                % compute the solution associated to the initial condition
                lowerThanKc = (IC(:,3) <= k_c);
                greaterThanKc = (IC(:,3) > k_c);
                
                %INITIAL==== <k_c ===========================================
                % first check those conditions with <= k_c
                inCharacterDomain = IC(:,1)+v_f*time <= position &...
                                    IC(:,2)+v_f*time >= position;
                
                % only applies to those <k_c
                activeCharacterDomain = lowerThanKc & inCharacterDomain;
                
                inFanDomain = IC(:,1)+v_f*time > position &...
                              IC(:,1)+w*time <= position;
                
                activeFanDomain = lowerThanKc & inFanDomain;
                
                % solution in characteristic domain
                tmp_M = min( IC_cum_num_veh(activeCharacterDomain,1) -...
                    IC(activeCharacterDomain,3).*(position - v_f*time - ...
                    IC(activeCharacterDomain,1)) );
                M(i) = self.minNonEmpty(M(i), tmp_M);
                % solution in fan domain
                tmp_M = min( IC_cum_num_veh(activeFanDomain,1) -...
                    k_c.*(position - v_f*time - ...
                    IC(activeFanDomain,1)  )  );
                M(i) = self.minNonEmpty(M(i), tmp_M);


                %INITIAL==== >k_c ===========================================
                % now check the initial condition with > k_c
                inCharacterDomain = IC(:,1)+w*time <= position &...
                                    IC(:,2)+w*time >= position; 
                activeCharacterDomain = greaterThanKc & inCharacterDomain;
                
                inFanDomain = IC(:,2)+v_f*time > position &...
                              IC(:,2)+w*time <= position; 
                activeFanDomain = greaterThanKc & inFanDomain;
                
                % solution in characteristic domain
                tmp_M = min( IC_cum_num_veh(activeCharacterDomain,1) -...
                    IC(activeCharacterDomain,3).*(position - w*time - ...
                    IC(activeCharacterDomain,1)) -k_m*time*w);
                M(i) = self.minNonEmpty(M(i), tmp_M);
                % solution in fan domain
                tmp_M = min( IC_cum_num_veh(activeFanDomain,1)+ IC_num_veh(activeFanDomain,1)...
                    - k_c.*(position - w*time - IC(activeFanDomain,2)) ...
                    - k_m*time*w);
                M(i) = self.minNonEmpty(M(i), tmp_M);

                %UPSTREAM==================================================
                % now check the upstream conditions, only apply those
                % boundary conditions that are not NaN
                notNaN = ~isnan(BC_us(:,3));
                
                inCharacterDomain = BC_us(:,1)+ position/v_f <= time &...
                                    BC_us(:,2)+ position/v_f >= time; 
                activeCharacterDomain = notNaN & inCharacterDomain;
                
                inFanDomain = BC_us(:,2)+ position/v_f < time;
                activeFanDomain = notNaN & inFanDomain;

                % solution in characteristic domain
                tmp_M = min( BC_cum_us_M(activeCharacterDomain,1) +...
                    BC_us(activeCharacterDomain,3).*(time - position/v_f - ...
                    BC_us(activeCharacterDomain,1)) );
                M(i) = self.minNonEmpty(M(i), tmp_M);
                % solution in fan domain
                tmp_M = min( BC_cum_us_M(activeFanDomain,1) +...
                    BC_us_num_veh(activeFanDomain,1) + ...
                    k_c*v_f.*(time - position/v_f - ...
                    BC_us(activeFanDomain,2)) );
                M(i) = self.minNonEmpty(M(i), tmp_M);

                %DOWNSTREAM==================================================
                % now check the upstream conditions
                notNaN = ~isnan(BC_ds(:,3));
                
                inCharacterDomain = BC_ds(:,1) + (position - len_link)/w <= time &...
                                    BC_ds(:,2)+ (position - len_link)/w >= time;
                activeCharacterDomain = notNaN & inCharacterDomain;
                
                inFanDomain = BC_ds(:,2)+ (position - len_link)/w < time;
                activeFanDomain = notNaN & inFanDomain;
                
                % solution in characteristic domain
                tmp_M = min( BC_cum_ds_M(activeCharacterDomain,1) +...
                    BC_ds(activeCharacterDomain,3).*(time -...
                    (position - len_link)/w - ...
                    BC_ds(activeCharacterDomain,1))...
                    - k_m*(position - len_link));
                M(i) = self.minNonEmpty(M(i), tmp_M);
                % solution in fan domain
                tmp_M = min( BC_cum_ds_M(activeFanDomain,1) +...
                    BC_ds_num_veh(activeFanDomain,1) + ...
                    k_c*v_f.*(time - (position-len_link)/v_f - ...
                    BC_ds(activeFanDomain,2)) );
                
                M(i) = self.minNonEmpty(M(i), tmp_M);


            end
        
        end



        %===============================================================
        function pt_found = searchShocksOnLine(self, link, t_interval, pos_interval, ...
                                                    M_interval, searchDepth, tol)
            % This function returns the intersection of shock waves with a line defined by the two points, not including those two points.
            %   - it can be used to find the intersection point with horizontal or vertical line segment.
            %   - View the vehicle id as a function of time and find the corner point.
            %   - If it is a (or an approximated < dt_tol) vertical line, then 
            %    view vehicle id as a function of the position, then find the corner point.
            %   - Once the time or the position of the corner point found, then find 
            %    the correspoing position or time on the line.
            % input:
            %       link: the link ID
            %       t_interval: 2x1 float, the time of two points
            %       pos_interval: 2x1 float, the position of two points in meters
            %       M_interval: 2x1 float, the absolute vehicle ID, to save computation time
            %                   set NaN if not sampled yet
            %       searchDepth: int, the maximal search depth of the recursive function
            %       tol: 2x1 float [dt_tol; dx_tol] the minimal value that we would further
            %            split and locate the intersection point
            % output:
            %       pt_found: struct, .time, a column vector
            %                         .position, a column vector, in meters
            %                         .vehicle_id, a column vector, absolute vehicle ID
            
            pt_found = [];

            % if search depth is 
            if searchDepth > self.max_search_depth
                return 
            end

            % First, compute the vheicle id at two bounds
            t_L = t_interval(1);
            t_R = t_interval(2);
            t_C = (t_L+t_R)/2;
            pos_L = pos_interval(1);
            pos_R = pos_interval(2);
            pos_C = (pos_L+pos_R)/2;
            % If value not computed yet
            if isnan(M_interval(1))
                M_L = self.absVehID(link, t_L, pos_L);
            else
                M_L = M_interval(1);
            end
            if isnan(M_interval(2))
                M_R = self.absVehID(link, t_R, pos_R);
            else
                M_R = M_interval(2);
            end
            M_C = self.absVehID(link, t_C, pos_C);
            
            % Second, check if they are on a straight line
            % if t_interval > dt_tol, then find intersection by time
            if abs(t_interval(2) - t_interval(1)) > tol(1)

                if self.onStraightLine([t_L; t_C; t_R], [M_L; M_C; M_R])

                    % double check if on a straight line. Some times, the
                    % segments are symmetric and those 3 points are on the same
                    % line, but the sending function is not a line
                    t_Q_1 = t_L + (t_R-t_L)/4;      % first quarter
                    t_Q_3 = t_L + 3*(t_R-t_L)/4;    % third quarter
                    pos_Q_1 = pos_L + (pos_R-pos_L)/4;
                    pos_Q_3 = pos_L + 3*(pos_R-pos_L)/4;
                    
                    M_Q_1 = self.absVehID(link, t_Q_1, pos_Q_1);
                    M_Q_3 = self.absVehID(link, t_Q_3, pos_Q_3);
                    
                    if self.onStraightLine([t_L; t_Q_1; t_C], [M_L; M_Q_1; M_C]) &&...
                       self.onStraightLine([t_C; t_Q_3; t_R], [M_C; M_Q_3; M_R])
                        % sampled 5 points on the same line, then very likely
                        % to be a straight line
                        return
                    end
                end

                % Third, if not on a straing line, then find the left and right slope
                % because the slope may be computed using time or position, we do not 
                % pass the last computed value in recursion.
                s_L = self.findFuncSlope(link, [t_L; t_C], [pos_L; pos_C], ...
                                         [M_L; M_C], 'left', 'time', tol);
                s_R = self.findFuncSlope(link, [t_C; t_R], [pos_C; pos_R], ...
                                         [M_C; M_R], 'right', 'time', tol);

                % compute the intersection time
                if abs(s_L - s_R) >= 1.0e-3
                    % make sure two slopes are not parallel
                    t_insct = (M_R - M_L + s_L*t_L - s_R*t_R)/(s_L - s_R);
                    M_insct = s_L*(t_insct - t_L) + M_L;
                    
                    if t_insct <= t_R && t_insct >= t_L
                        % if the intersection is in the time step, then
                        % validate, otherwise, move on to split
                        pos_insct = (pos_R-pos_L)*(t_insct-t_L)/(t_R-t_L) + pos_L;
                        M_validate = self.absVehID(link, t_insct, pos_insct);
                        if abs(M_insct - M_validate) <= 1e-6
                            % found the intersection
                            pt_found = [pt_found; t_insct, pos_insct, M_insct];
                            return
                        end
                    end
                end

                % Forth, recursion, split the interval if there are more than 
                % one corners based on: 
                % 1. slopes are parallel; 2. intersection outside of
                % the interval; 3. cornerpoint not valid in sample
                pt_found_left = self.searchShocksOnLine(link, [t_L; t_C], [pos_L; pos_C],...
                                    [M_L; M_C], searchDepth+1, tol);
                pt_found_right = self.searchShocksOnLine(link, [t_C; t_R], [pos_C; pos_R],...
                                    [M_C; M_R], searchDepth+1, tol);
                pt_found = [pt_found; pt_found_left; pt_found_right];


            % otherwise find by position
            elseif abs(pos_interval(2) - pos_interval(1)) > tol(2)
                
                if self.onStraightLine([pos_L; pos_C; pos_R], [M_L; M_C; M_R])

                    % double check if on a straight line. Some times, the
                    % segments are symmetric and those 3 points are on the same
                    % line, but the sending function is not a line
                    t_Q_1 = t_L + (t_R-t_L)/4;      % first quarter
                    t_Q_3 = t_L + 3*(t_R-t_L)/4;    % third quarter
                    pos_Q_1 = pos_L + (pos_R-pos_L)/4;
                    pos_Q_3 = pos_L + 3*(pos_R-pos_L)/4;
                    
                    M_Q_1 = self.absVehID(link, t_Q_1, pos_Q_1);
                    M_Q_3 = self.absVehID(link, t_Q_3, pos_Q_3);
                    
                    if self.onStraightLine([pos_L; pos_Q_1; pos_C], [M_L; M_Q_1; M_C]) &&...
                       self.onStraightLine([pos_C; pos_Q_3; pos_R], [M_C; M_Q_3; M_R])
                        % sampled 5 points on the same line, then very likely
                        % to be a straight line, return []
                        return
                    end
                end

                % Third, if not on a straing line, then find the left and right slope
                % because the slope may be computed using time or position, we do not 
                % pass the last computed value in recursion.
                s_L = self.findFuncSlope(link, [t_L; t_C], [pos_L; pos_C], ...
                                         [M_L; M_C], 'left', 'position', tol);
                s_R = self.findFuncSlope(link, [t_C; t_R], [pos_C; pos_R], ...
                                         [M_C; M_R], 'right', 'position', tol);

                % compute the intersection time
                if abs(s_L - s_R) >= 1.0e-3
                    % make sure two slopes are not parallel
                    pos_insct = (M_R - M_L + s_L*pos_L - s_R*pos_R)/(s_L - s_R);
                    M_insct = s_L*(pos_insct - pos_L) + M_L;
                    
                    if pos_insct <= pos_R && pos_insct >= pos_L
                        % if the intersection is in the time step, then
                        % validate, otherwise, move on to split
                        t_insct = (t_R-t_L)*(pos_insct-pos_L)/(pos_R-pos_L) + t_L;
                        M_validate = self.absVehID(link, t_insct, pos_insct);
                        if abs(M_insct - M_validate) <= 1e-6
                            % found the intersection
                            pt_found = [pt_found; t_insct, pos_insct, M_insct];
                            return
                        end
                    end
                end

                % Forth, recursion, split the interval if there are more than 
                % one corners based on: 
                % 1. slopes are parallel; 2. intersection outside of
                % the interval; 3. cornerpoint not valid in sample
                pt_found_left = self.searchShocksOnLine(link, [t_L; t_C], [pos_L; pos_C],...
                                    [M_L; M_C], searchDepth+1, tol);
                pt_found_right = self.searchShocksOnLine(link, [t_C; t_R], [pos_C; pos_R],...
                                    [M_C; M_R], searchDepth+1, tol);
                pt_found = [pt_found; pt_found_left; pt_found_right];


            % if the area is too small to further split for more refined search 
            % then just return the center point
            else
                
                pt_found = [sum(t_interval)/2, sum(pos_interval)/2,...
                            self.absVehID(link, sum(t_interval)/2, sum(pos_interval)/2)];
                return

            end

            

        end


        %===============================================================
        function density = extractDensity(self, time)
            % This function returns the density estimates at a given time.
            %   - It only computes the vehicle IDs at that time, hence no need 
            %    to compute the entire density estimation diagram, faster.
            %   - it searches shockwave intersections at that time, and the 
            %    traffic density will be aggregated and separated by the 
            %    intersectoin points. Hence less initial conditions, faster.
            % input: 
            %       time: float, the time that we would like to extract density
            % output:
            %       density: struct, .(linkStr).X_grid_cum
            %                        .(linkStr).IC
            
            % specify the search depth and tolerance
            search_depth = 0;
            tol = [1e-3; 1e-2];
            
            % check each link
            for link = self.net.link_labels'

                linkStr = sprintf('link_%d', link);

                len_link = self.net.network_hwy.(linkStr).para_postkm*1000;

                pt_found = self.searchShocksOnLine(link, [time; time], ...
                    [0; len_link], [NaN; NaN], search_depth, tol);

                % get the time, position, and vehicle id at two bounds
                pt_us = [time, 0, self.absVehID(link, time, 0)];
                pt_ds = [time, len_link, self.absVehID(link, time, len_link)];

                % all the grid points
                pt_all = [pt_us; pt_found; pt_ds];
                pt_grid.position = pt_all(:,2);
                pt_grid.vehicle_id = pt_all(:,3);

                % compute the density
                dx = pt_grid.position(2:end) - pt_grid.position(1:end-1);
                dM = pt_grid.vehicle_id(1:end-1) - pt_grid.vehicle_id(2:end);
                rho = dM./dx;

                % normalize to kc
                density.(linkStr).IC = rho;
                density.(linkStr).X_grid_cum = pt_grid.position;

            end


        end

        

        %===============================================================
        function slope = findFuncSlope(self, link, t_interval, pos_interval, ...
                                       M_interval, slope_side, wrt, tol)
            % This funciton returns the slope of the vehicle id with respect to time or position.
            % This is different from findBoundaryFuncSlope function since findBoundaryFuncSlope uses 
            % sendingFuncAtLink which removes the cloese boundary condition, while this function uses 
            % absVehID which get the vehicle id at a certain points
            % input:
            %       link: the link id
            %       t_interval: 2x1 float vector
            %       pos_interval: 2x1 float vector define the line
            %       M_interval: 2x1 float vector, the absolute vehicle id
            %       slope_side: 'left', the first point side; 'right', the second point side
            %       wrt: 'time', 'position', with respect to
            %       tol: [dt_tol; dx_tol], the tolerance for stopping recursion.
            
            if strcmp(wrt, 'time')

                % if the time interval is too small, then stop computing slope since 
                % it may introduce numerical error
                if t_interval(2)-t_interval(1) <= tol(1)
                    warning('WARNING: The slope can not be computed\n')
                    slope = NaN;
                    return
                end

                % First, sample boundary points values and the center point
                t_L = t_interval(1);
                t_R = t_interval(2);
                t_C = (t_L+t_R)/2;
                pos_L = pos_interval(1);
                pos_R = pos_interval(2);
                pos_C = (pos_L+pos_R)/2;
                % If value not computed yet
                if isnan(M_interval(1))
                    M_L = self.absVehID(link, t_L, pos_L);
                else
                    M_L = M_interval(1);
                end
                if isnan(M_interval(2))
                    M_R = self.absVehID(link, t_R, pos_R);
                else
                    M_R = M_interval(2);
                end
                
                M_C = self.absVehID(link, t_C, pos_C);

                if self.onStraightLine([t_L; t_C; t_R],[M_L; M_C; M_R])
                    % if true, return the slope
                    slope = (M_R-M_L)/(t_R-t_L);
                else
                    if strcmp(slope_side, 'left')
                        % keep the left half to find the slope
                        slope = self.findFuncSlope(link, [t_L; t_C], [pos_L; pos_C],...
                                         [M_L; M_C], slope_side, wrt, tol);
                    elseif strcmp(slope_side, 'right')
                        % keep the right half to find the slope
                        slope = self.findFuncSlope(link, [t_C; t_R], [pos_C; pos_R],...
                                         [M_C; M_R], slope_side, wrt, tol);
                    end
                end                           



            elseif strcmp(wrt, 'position')
                
                if pos_interval(2)-pos_interval(1) <= tol(2)
                    warning('WARNING: The slope can not be computed\n')
                    slope = NaN;
                    return
                end

                % First, sample boundary points values and the center point
                t_L = t_interval(1);
                t_R = t_interval(2);
                t_C = (t_L+t_R)/2;
                pos_L = pos_interval(1);
                pos_R = pos_interval(2);
                pos_C = (pos_L+pos_R)/2;
                % If value not computed yet
                if isnan(M_interval(1))
                    M_L = self.absVehID(link, t_L, pos_L);
                else
                    M_L = M_interval(1);
                end
                if isnan(M_interval(2))
                    M_R = self.absVehID(link, t_R, pos_R);
                else
                    M_R = M_interval(2);
                end
                
                M_C = self.absVehID(link, t_C, pos_C);

                if self.onStraightLine([pos_L; pos_C; pos_R],[M_L; M_C; M_R])
                    % if true, return the slope
                    slope = (M_R-M_L)/(pos_R-pos_L);
                else
                    if strcmp(slope_side, 'left')
                        % keep the left half to find the slope
                        slope = self.findFuncSlope(link, [t_L; t_C], [pos_L; pos_C],...
                                         [M_L; M_C], slope_side, wrt, tol);
                    elseif strcmp(slope_side, 'right')
                        % keep the right half to find the slope
                        slope = self.findFuncSlope(link, [t_C; t_R], [pos_C; pos_R],...
                                         [M_C; M_R], slope_side, wrt, tol);
                    end
                end      



            end

        end


        %===============================================================
        function [k_trans] = mapping(~, k_ori, original_interval, target_interval)
            % This function maps the values from original interval to a target interval
            % input:
            %       k_ori: the original values to be mapped
            %       original_interval: nx2 matrix. Each row is an interval to be mapped
            %           to the corresponding interval in target_interval.
            %       target_interval: nx2 matrix. Each row is the target interval.
            % output:
            %       k_trans: the mapped values which has the same dimension as k_ori
            % Example:
            % k_trans = mapping(k,[0 k_c; k_c k_m],[0 0.5*k_m; 0.5*k_m k_m]);
            % mapping (0~k_c)(k_c~k_,m) ==> (0~0.5*k_m)(0.5*k_m~k_m)
            
            k_trans = 0.0*k_ori;
            for i = 1:size(original_interval,1)
                isInDomain = (k_ori>=original_interval(i,1) & k_ori<=original_interval(i,2)+10e-6);
                k_trans(isInDomain) = target_interval(i,1) + ...
                    (k_ori(isInDomain)-original_interval(i,1))*(target_interval(i,2)-target_interval(i,1))...
                    /(original_interval(i,2)-original_interval(i,1));
            end
            
            
        end
        
        
        
        %===============================================================
        function v = minNonEmpty(~, v1, v2)
            % The function compares the min of two values which may be NaN. 
            % If one value is NaN, it will return the other value. If both NaN, return NaN.
            
            if isempty(v1) && ~isempty(v2)
                v = min(v2);
            end
            
            if isempty(v2) && ~isempty(v1)
                v = min(v1);
            end
            
            if isempty(v1) && isempty(v2)
                v = NaN;
            end
            
            if ~isempty(v1) && ~isempty(v2)
                v = min( min(v1),min(v2) );
            end
            
        end

        
        
        %===============================================================
        function [uniq] = unique_tol(~, array, tol)
            % This function has the same functionality as unique() but considers two values within tol as identical. 
            % input:
            %       array: the array that we would like to find the unique values.
            %       tol: the tolerance between two unique values.
            % output: 
            %       uniq: the unique values in the array
            
            uniq = unique(array);
            
            len = length(uniq);
            index = (uniq(2:len)-uniq(1:len-1))<tol;
            
            % remove the value
            if sum(index)~=0
                uniq(index) = [];
            end
            
        end
        
        
        
    end     %end methods


end

























