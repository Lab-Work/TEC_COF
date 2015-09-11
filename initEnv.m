% initEnv
% initialize the simulation environment, including parameters, index


%% Path
% clear global path path_data path_profiles path_result;
% global path path_data path_result;
% path = [pwd '/'];   
% path_data = [path 'data/'];    % place to input data
% path_result = [path 'result/'];
% if ~exist(path_result, 'dir')
%     mkdir(path_result);
% end

%% default_para
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
default_para.e_max = 0.0;









