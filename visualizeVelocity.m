%% set up Grid, Initial Data
%how many grid points we want in all directions
gpoints = 45;
N = gpoints*ones(4,1);         % Number of grid points per dimension, default to 41

%bounds on all possible grid points
grid_max = -[-1.8091-pi/15, -0.3289-pi/15, 0.0013-pi/15, -4.9904-pi/15];
grid_min = -[0.0247+pi/15, 4.4269+pi/15, 2.4655+pi/15, 1.8952+pi/15];

%g = createGrid(grid_min, grid_max, N);

%set what counts as "standing"
%data = shapeRectangleByCorners(g, lowerlimit, upperlimit)
%data = shapeRectangleByCorners(g, -[pi/80;.1; pi/15;.1], -[-pi/15;-.1; -pi/80;-.1]); 

%% parameters
height = 1.72;
mass = 80;
M1 = 2*(0.1416*mass);       % mass of thighs 
M2 = (.0694 +.4346)*mass;    % mass of head-arms-trunk
L0 = .25*height;          % length of segment (shank)
L1 = .26*height;          % length of segment .4(thigh)
R1 = .43*L1;          % position of COM along segment (thigh)
R2 = .6*.4*height;
L2 = .3*height;

T1Max = 107;
T1Min = -T1Max;
T2Max = 87;
T2Min = -60;
TAMax = 68;
TAMin = -50;

alpha = 1;
T1Max = alpha*T1Max;
T1Min = alpha*T1Min;
T2Max = alpha*T2Max;
T2Min = alpha*T2Min;
TAMax = alpha*TAMax;
TAMin = alpha*TAMin;



%% put all this in schemeData
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

%Use testAnkle code to define allowable torques, save to schemeData
schemeData = testAnkle(schemeData);

%% finding qualified points within constraints
A = schemeData.tau1test{1}; %pull out grid of points that have feasible points to be tested

b = 1;
%for each of the 8 possible points
for m = 2:8
  A = min(A,schemeData.tau1test{m}); %find all the states the have allowable torques
  %by taking the min across all of these points we're throwing out any
  %points that only have "NaN" as possible torques
end

% we care only about the points on this grid that have feasible torques.
% we're going to do some awkward manipulation to make sure we only have
% those:
A=(A.*0)+1;
data = data.*A;

% Anan = isnan(A); %find points that are not allowed on 2D grid
% Anew = Anan<1;
% A = Anew.+

% 
% gtest1 = g2D.xs{1}.*Anew; %removed not-allowed points from grid
% gtest2 = g2D.xs{2}.*Anew;
% 
% gtest1_rem = g2D.xs{1}.*Anan; %keep not-allowed points from grid
% gtest2_rem = g2D.xs{2}.*Anan;

%slice = -pi/2;

B=data.*(data>-100000);
[g2D,data2D]=proj(g,data,[1 0 1 0],'min');

h = visSetIm(g2D, data2D); %, color, level, extraArgs);
