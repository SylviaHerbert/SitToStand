function [data, g, tau, schemeData] = Double_KneeHips_Backward(gpoints, accuracy, tMax, alpha,schemeData)
%% Grid
grid_min = [-0.29, -6.31, -2.67, -7.57];
grid_max = [1.89, 4.51, 0.15, 8.91];
N = gpoints*ones(4,1);  
g = createGrid(grid_min, grid_max, N);


%Want to drop # of test points from 8 to 6? use trim.
trim = 1;
%% target set

max_v = (pi/8);       % allowing for some sway
standing_min = [-pi/15, -max_v, -pi/15, -max_v];
standing_max = [pi/15, max_v, 0.15, max_v];

data = shapeRectangleByCorners(g, standing_min, standing_max);

%% time vector
t0 = 0;

if nargin < 3
tMax = 1;
end

dt = 0.01;
tau = t0:dt:tMax;
% If intermediate results are not needed, use tau = [t0 tMax];

%% problem parameters
% T1Max = 1;
% T2Max = 1;
height = 1.72;
mass = 62;
M1 = 2*(0.1416*mass);       % mass of thighs 
M2 = (.0694 +.4346)*mass;    % mass of head-arms-trunk
L0 = .25*height;          % length of segment (shank)
L1 = .26*height;          % length of segment .4(thigh)
R1 = .43*L1;          % position of COM along segment (thigh)
R2 = .6*.4*height;          % position of COM along segment (head-arms-trunk)


  T1Max = 107;
T1Min = -T1Max;
T2Max = 87;
T2Min = -60;
TAMax = 68;
TAMin = -50;
if nargin < 4
alpha = 1;
end


T1Max = alpha*T1Max;
T1Min = alpha*T1Min;
T2Max = alpha*T2Max;
T2Min = alpha*T2Min;
TAMax = alpha*TAMax;
TAMin = alpha*TAMin;
%% Pack problem parameters
if nargin <5
schemeData.grid = g; % Grid MUST be specified!
schemeData.T1Max = T1Max;
schemeData.T1Min = T1Min;
schemeData.T2Max = T2Max;
schemeData.T2Min = T2Min;
schemeData.TAMax = TAMax;
schemeData.TAMin = TAMin;
schemeData.R1 = R1;
schemeData.R2 = R2;
schemeData.M1 = M1;
schemeData.M2 = M2;
schemeData.L1 = L1;
schemeData.L0 = L0;
schemeData.height = height;
schemeData.accuracy = accuracy; %default is medium
schemeData = testAnkle(schemeData,trim); %add in ankle constraints
end

% ----- System dynamics are specified here -----
schemeData.hamFunc = @pendulum4DHam;
schemeData.partialFunc = @pendulum4Dpartial;
%% solve
%extraArgs.visualize = true;
%extraArgs.save_filename = 'Double_KneeHips_Backward_41_highAc';
%extraArgs.saveFrequency = 100;
%extraArgs.stopInit = [pi/2 0 -pi/2 0];
extraArgs.stopConverge = 1;
[data, tau] = HJIPDE_solve( ...
  data, tau, schemeData, 'zero', extraArgs);
%data = min(data,[],5);
end
