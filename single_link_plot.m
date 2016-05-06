% Single link solver and plot

time_grid = [0, 20, 40, 60, 80, 100];
q_in = [0.8,1.2,0.2,1.4,1.0];
q_out = [0.3,0.6,1.4,0.6,0.8];

ini_seg = [0,250,500,750,1000];
p_ini = [0.01,0.1,0.2,0.03];

% in m, s
vf = 24.6;
w = -4;
km = 0.4053;
kc = 0.0567;

us_position = 0;
ds_position = 1000;

% 62.5 x 25 matrix
dx_res = 40;
dt_res = 1.6;

% high resolution
dx_res = 1;
dt_res = 0.1;

% ==============================
fd = LH_Tfd(vf,w,km);

pbEnv = LH_general(fd,us_position,ds_position);


pbEnv.setIniDens(ini_seg,p_ini);
pbEnv.setUsFlows(time_grid,q_in);
pbEnv.setDsFlows(time_grid,q_out);


x_mesh_m = 0:dx_res:ds_position;
t_mesh_s = 0:dt_res:100;

xValues = ones(size(t_mesh_s'))*(x_mesh_m);
tValues = t_mesh_s' * ones(size(x_mesh_m));

result = pbEnv.explSol(tValues,xValues);
                
N = result{1};
activeComp = result{2};
k = pbEnv.density(tValues,xValues,activeComp);

% plot
%===========================================
%Transformation for better color presentation
%kc=>0.5km, km=>km
k_c_tmp = kc;
k_m_tmp = km;
k_trans = mapping(k, [0 k_c_tmp; k_c_tmp k_m_tmp],...
    [0 0.5*k_m_tmp; 0.5*k_m_tmp k_m_tmp]);

%===========================================
scrsz = get(0,'ScreenSize');
figure('Position',[1 1 scrsz(3) scrsz(4)]);
title(sprintf('Link'),'fontsize',24);

colormap jet

[~, ~] = LH_plot2D(t_mesh_s, x_mesh_m, N, k_trans, fd);


scrsz = get(0,'ScreenSize');
figure('Position',[1 1 scrsz(3) scrsz(4)]);
image(flipud(k_trans'),'CDataMapping','scaled')
colormap jet


















