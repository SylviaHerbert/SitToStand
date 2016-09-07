function [g2Dforward, data2Dforward, data0] = Single_Hips_Forwards(T2Max)
%% Grid
grid_min = [0.0013-pi/15, -4.9904-pi/15]; % Lower corner of computation domain
grid_max = [2.4655+pi/15, 1.8952+pi/15];    % Upper corner of computation domain
N = [41; 41];         % Number of grid points per dimension
%pdDims = [1,3];               % 1st, 3rd dimension is periodic
g2Dforward = createGrid(grid_min, grid_max, N);
% Use "g = createGrid(grid_min, grid_max, N);" if there are no periodic
% state space dimensions

%% target set

data0 = shapeRectangleByCorners(g2Dforward, [-pi/80;-.1],[pi/15;.1]);

% also try shapeRectangleByCorners, shapeSphere, etc.

%% time vector
t0 = 0;
tMax = 1;
dt = 0.025;
tau = t0:dt:tMax;
% If intermediate results are not needed, use tau = [t0 tMax];

%% problem parameters
%T2Max = 1;
R2 = .222;
M2 = 27.3;

%% Pack problem parameters
schemeData.grid = g2Dforward; % Grid MUST be specified!
schemeData.T2Max = T2Max;
schemeData.R2 = R2;
schemeData.M2 = M2;

% ----- System dynamics are specified here -----
schemeData.hamFunc = @pendulum2DHam;
schemeData.partialFunc = @pendulum2DPartial;
%% solve
extraArgs.visualize = true;
[data2Dforward, tau2D, extraOuts2D] = HJIPDE_solve( ...
  data0, tau, schemeData, 'zero',extraArgs);
end