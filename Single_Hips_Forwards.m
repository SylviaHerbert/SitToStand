%% Grid
grid_min = [.0013-.5; -4.99-.5]; % Lower corner of computation domain
grid_max = [2.4655+.5; 1.8952+.5];    % Upper corner of computation domain
N = [41; 41];         % Number of grid points per dimension
%pdDims = [1,3];               % 1st, 3rd dimension is periodic
g2D = createGrid(grid_min, grid_max, N);
% Use "g = createGrid(grid_min, grid_max, N);" if there are no periodic
% state space dimensions

%% target set

data0 = shapeRectangleByCorners(g2D, [-pi/80;-.1],[pi/15;.1]);

% also try shapeRectangleByCorners, shapeSphere, etc.

%% time vector
t0 = 0;
tMax = 2;
dt = 0.025;
tau = t0:dt:tMax;
% If intermediate results are not needed, use tau = [t0 tMax];

%% problem parameters
T2Max = 40;
R2 = .222;
M2 = 27.3;

%% Pack problem parameters
schemeData.grid = g2D; % Grid MUST be specified!
schemeData.T2Max = T2Max;
schemeData.R2 = R2;
schemeData.M2 = M2;

% ----- System dynamics are specified here -----
schemeData.hamFunc = @pendulum2DHam;
schemeData.partialFunc = @pendulum2DPartial;
%% solve
extraArgs.visualize = true;
[data2D, tau2D, extraOuts2D] = HJIPDE_solve( ...
  data0, tau, schemeData, 'zero',extraArgs);