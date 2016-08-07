
%% Grid
grid_min = [-1.8091-.5; -3.289-.5; .0013-.5; -4.99-.5]; % Lower corner of computation domain
grid_max = [.0247+.5; 4.4269+.5; 2.4655+.5; 1.8952+.5];    % Upper corner of computation domain
N = [41; 41; 41; 41];         % Number of grid points per dimension
%pdDims = [1,3];               % 1st, 3rd dimension is periodic
g = createGrid(grid_min, grid_max, N);
% Use "g = createGrid(grid_min, grid_max, N);" if there are no periodic
% state space dimensions

%% target set

data0 = shapeRectangleByCorners(g, [-pi/15;-.1; -pi/80;-.1],[pi/80;.1; pi/15;.1]);

% also try shapeRectangleByCorners, shapeSphere, etc.

%% time vector
t0 = 0;
tMax = 2;
dt = 0.025;
tau = t0:dt:tMax;
% If intermediate results are not needed, use tau = [t0 tMax];

%% problem parameters
T1Max = 40;
T2Max = 40;
R1 = .276;
R2 = .222;
M1 = 12.4;
M2 = 27.3;
L1 = .438;

%% Pack problem parameters
schemeData.grid = g; % Grid MUST be specified!
schemeData.T1Max = T1Max;
schemeData.T2Max = T2Max;
schemeData.R1 = R1;
schemeData.R2 = R2;
schemeData.M1 = M1;
schemeData.M2 = M2;
schemeData.L1 = L1;

% ----- System dynamics are specified here -----
schemeData.hamFunc = @pendulum4DHam;
schemeData.partialFunc = @pendulum4Dpartial;
%% solve
%extraArgs.visualize = true;
[data, tau, extraOuts] = HJIPDE_solve( ...
  data0, tau, schemeData, 'zero');