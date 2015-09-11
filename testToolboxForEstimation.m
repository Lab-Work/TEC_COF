% This script is a simple forward simulation for debugging the COF toolbox. 
% This is an example of three links, which is a merge
% Yanning, Sep 04, 2015

% This updated script is to test the modified toolbox in COF_control
% folder.

tic
clearvars -except dbg

profile on
%===============================================================
% configuration paramters
runCtrl = 'no_control'; % This parameter can be later used for optimal control 
entropyTolerance = 0.2; % number of vehicles that we would consider as entropic solution

sim_steps = 5;
step_length = 30;   % seconds
T_BC_data = ones(sim_steps,1)*step_length;
T_BC_cum = [0; cumsum(T_BC_data)];


start_time = 0;
end_time = sum(T_BC_data);

%===============================================================
% This is the time discretization grid. Will be iteratively updated to
% eliminate the discretization error.
T_junc= T_BC_data;

% initialize environment, parameters, index.
initEnv;    

%===============================================================
if strcmp(runCtrl,'no_control') == 1
    % here we define an example
    % Maximize the out flow of a single link

    % test length
    len_link1 = 1;    %km
    len_link2 = 1;
    len_link3 = 1;
    
    % Initial traffic density is initialized as segments with equal length
    % Normalized to rho_c
%     rho_tmp = randn(sim_steps, 1) + 0.5;
%     rho_tmp(rho_tmp <= 0.2) = 0.2;
%     rho_tmp(rho_tmp >= 1.5) = 1.5;
%     Ini_1 = rho_tmp;
    Ini_1.IC = [2, 2, 1, 0.2, 0.2]';

%     rho_tmp = randn(sim_steps, 1) + 0.5;
%     rho_tmp(rho_tmp <= 0.2) = 0.2;
%     rho_tmp(rho_tmp >= 1.5) = 1.5;
%     Ini_2 = rho_tmp;
    Ini_2.IC = [1, 0.5, 1, 0.5, 1]';
    Ini_3.IC = zeros(sim_steps, 1);
    
    % Boundary condition at the two entrances and the one exit
    % Normalized to q_max
    % generate a random upstream flow
%     q_tmp = randn(sim_steps, 1) + 1;
%     q_tmp(q_tmp <= 0.5) = 0.5;
%     q_tmp(q_tmp >= 0.9) = 0.9;
%     q1_us_data = q_tmp;
    q1_us_data = [0.5, 0.9, 0.5, 0.9, 0.5]';
    q2_us_data = [0.9, 0.5, 0.9, 0.5, 0.9]';
    
%     q_tmp = randn(sim_steps, 1) + 1;
%     q_tmp(q_tmp <= 0.5) = 0.5;
%     q_tmp(q_tmp >= 0.9) = 0.9;
%     q2_ds_data = q_tmp;
    
end

% represent the boundary condition with the time grid information
q1_us_TF = [T_BC_cum(1:sim_steps), T_BC_cum(2:sim_steps+1), q1_us_data];
q2_us_TF = [T_BC_cum(1:sim_steps), T_BC_cum(2:sim_steps+1), q2_us_data];

% This is the while loop for iteratively regridding and eventually getting
% the entropy condition
% Flags for iteratively updating time discretization
getEntropy = false;
loopCounter = 0;
while getEntropy == false && loopCounter <=50
    
    loopCounter = loopCounter+1;
    
    % Keep track of the objective funciton values, the solution, and the
    % split ratio
    f = [];
    Sol = [];
    
    %===============================================================
    % Define the network
    % initLink(obj, link, para, num_lanes, lengthKM, linkType)
    net = initNetwork;
    net.addLink(1, default_para, 1, len_link1, 'freeway');
    net.addLink(2, default_para, 1, len_link2, 'freeway');
    net.addLink(3, default_para, 1, len_link3, 'freeway');
    
    % function addJunc(self, junc, inlabel, outlabel, type_junc, ratio, T)
    net.addJunc(1, [1, 2]', 3, 'merge', [1; 3], T_junc);
    
    %===============================================================                         
    % set initial conditoins
    net.setInitialCon(1,Ini_1);
    net.setInitialCon(2,Ini_2);
    net.setInitialCon(3,Ini_3);
    
    % setBoundaryCon(obj, link, q_in, q_out, T_in, T_out)
    net.setBoundaryCon(1, q1_us_data, [], T_BC_data, T_junc);
    net.setBoundaryCon(2, q2_us_data, [], T_BC_data, T_junc);
    net.setBoundaryCon(3, [], [], T_junc, T_BC_data);
    %===============================================================
    % define and solve optimization program
    LP = optProgram;
    LP.setConfig(start_time, end_time);
    LP.setConstraints(net, default_para.e_max);
    
    %===============================================================
    % Add objective functions
    LP.addEntropy(net, 1);
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
    Mos = postSolution(x, net, LP.dv_index, LP.end_time, 2, 2);
    Mos. estimateState();
    % Mos.plotLinks('all')
    Mos.plotJuncs('all')
    
    [getEntropy, steps] = Mos.checkEntropy(entropyTolerance);
    
    if getEntropy == false
        hold on
        T_junc_cum = Mos.updateTimeDiscretization(steps);
        % updated T_junc
        % a quick hack, only works for this particular one junc scenario
        T_junc = T_junc_cum.junc_1(2:end)-T_junc_cum.junc_1(1:end-1);
    end

end



toc


% profile viewer
% p = profile('info');
% profsave(p,'./result/matlab_profiler_results');





