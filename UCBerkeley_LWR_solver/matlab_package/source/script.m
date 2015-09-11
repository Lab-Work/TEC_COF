% In this example, all units are SI: meters, seconds.
% The algorithm itself is unit-agnostic, and only needs the unit used as
% input to be consistent.


% % This creates a triangular fundamental diagram
% % Free flow speed 30 m/s
% % Congestion wave speed -5 m/s
% % Maximum Density 0.1 veh / m
% % Uncomment line below if you want to create a triangulat fundamental
% % diagram
fd = LH_Tfd(30,-5,0.1);


%This creates a Greenshields fundamental diagram
% Free flow speed 30 m/s
% Maximum Density 0.1 veh / m
%fd = LH_Greenshields(30,.1);

%This creates a spatial domain between 0 and 1000, where the fundamental
%diagram is fd
%between 0 and 1000.
pbEnv = LH_general(fd,0,1000);

%Sample input arrays for initial densities.
%Initial densities: 80 veh/km at 0<x<200m, 10 veh/km at 200<x<500m,...
%pbEnv.setIniDens([0 500 1000], [80E-3 10E-3]);
%pbEnv.setIniDens([400 600], [80E-3]);

%Sample arrays for upstream flows.
%Upstream Flows: .4 veh/s at 0<t<20s, .01 veh/s at 20<t<40s...
%pbEnv.setUsFlows([0 20 40 50], [.4 .01 .2]);
% pbEnv.setUsFlows([0 35], [.2]);

%Sample arrays for downstream flows.
%Downstream flows: .3 veh/s at 0<t<30s, 0 veh/s at 30<t<35s...
%pbEnv.setDsFlows([0 30 35 50], [.3 0 .1]);
%pbEnv.setDsFlows([25 48], [.3]);
%pbEnv.setDsFlows([0 30 40], [.3 .1]);

%pbEnv.addIntCond(20,30,400,500,0.1);

pbEnv.addFirstIntCond(20,30,400,500,0.1);


%Here we create two matrices tValues and xValues which store space and time
%information for every point at which we want to compute the solution.
X = 1000;           % Maximal x in the computational domain
nx = 500;           % Number x grid points
T = 50;             % Maximal time for the computation
nt = 500;           % Number of t grid points
dx=X/nx;            % Space step
dt=T/nt;            % Time step
xScale = 0:dx:X;    % Create vector array for spatial domain
tScale = 0:dt:T;    % Create vector array for temporal domain

% Explication matrices pour donner l information (3 lignes max)
xValues = ones(size(tScale'))*(xScale);
tValues = tScale' * ones(size(xScale));

% tic/toc or [cputime] can be used to time the computation
tic
result = pbEnv.explSol(tValues,xValues);

% result{1} is the Moskowitz function 
N = result{1};

% result{2} is the active component matrix (needed for computation of density)
activeComp = result{2};

k = pbEnv.density(tValues,xValues,activeComp);%matrix of solution densities


LH_plot3D(tScale, xScale, N, k, fd)
figure
toc
LH_plot2D(tScale, xScale, N, k, fd);
