function [data, g, tau, schemeData, runtime] = Double_KneeHips_Backward(gpoints, accuracy, tMax, alpha, constraint, schemeData)

%% Grid
grid_min = [-0.29, -6.31, -2.67, -7.57];
grid_max = [1.89, 4.51, 0.15, 8.91];
gApoints = gpoints;
gVpoints = gApoints/2;
N = [gApoints gVpoints gApoints gVpoints]; 
g = createGrid(grid_min, grid_max, N);

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

dt = .01;
tau = t0:dt:tMax;

if nargin <4
  constraint = 0;
end
%% problem parameters

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

grav = 9.81;
%% Pack problem parameters
if nargin <6
  schemeData.dynSys= STS4D([0 0 0 0], R1, R2, M1, M2, L1, L0, grav,...
    T1Max, T1Min, T2Max, T2Min, TAMax, TAMin,[1:4]);
  schemeData.grid = g; % Grid MUST be specified!
  schemeData.accuracy = accuracy; %default is medium
  
  if constraint == 1
    %Want to drop # of test points from 8 to 6? use trim.
    trim = 1;
    [tau1s,tau2s, data] = testAnkle(schemeData, data, trim); %add in ankle constraints
    schemeData.dynSys.tau1Test = tau1s;
    schemeData.dynSys.tau2Test = tau2s;
  end

end

% ----- System dynamics are specified here -----
schemeData.hamFunc = @genericHam; %@pendulum4DHam;
schemeData.partialFunc = @genericPartial;%@pendulum4Dpartial;
%% solve
%extraArgs.visualize = true;
extraArgs.save_filename = (['Double_KneeHips_Backward_g', num2str(gpoints), ...
  '_', num2str(accuracy), '_t', num2str(tMax), '_alpha', num2str(alpha)]);
extraArgs.saveFrequency = 100;
extraArgs.stopInit = [pi/2 0 -pi/2 0];
extraArgs.stopConverge = 1;
extraArgs.keepLast = 1;
tic
[data, tau] = HJIPDE_solve( ...
  data, tau, schemeData, 'zero', extraArgs);
runtime = toc;
end
