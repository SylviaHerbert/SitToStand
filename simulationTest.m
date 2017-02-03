%say whether you want to see the allowed states given your torques, or just
%the standing conditions
clear; close all;
view = 'allowed_states';
% slice = [0 0];
 slice = 'min';


%% set up Grid, Initial Data
%how many grid points we want in all directions
gpoints = 45;
N = gpoints*ones(4,1);         % Number of grid points per dimension, default to 41

%bounds on all possible grid points
grid_max = -[-1.8091-pi/15, -0.3289-pi/15, 0.0013-pi/15, -4.9904-pi/15];
grid_min = -[0.0247+pi/15, 4.4269+pi/15, 2.4655+pi/15, 1.8952+pi/15];

g = createGrid(grid_min, grid_max, N);

%set what counts as "standing"
%data = shapeRectangleByCorners(g, lowerlimit, upperlimit)
data = shapeRectangleByCorners(g, -[pi/80;.1; pi/15;.1], -[-pi/15;-.1; -pi/80;-.1]); 

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

%% Project data through to 2D so we can view it

if strcmp(view,'allowed_states')
%%%NOTE!%%%
% if you want to see 2D through all possible velocities, use 'min'
% if you want to see 2D at a specific velocity, use [vKnee vHip], where 
% vKnee vHip are the velocities you want to use

[g2D, A2D]=proj(g,A,[0 1 0 1],slice); %project through to 2D
[g1D, A1D] = proj(g2D,A2D,[1 0],-pi/2);

% we care only about the points on this 2D grid that have feasible torques.
% we're going to do some awkward manipulation to make sure we only have
% those:
Anan = isnan(A2D); %find points that are not allowed on 2D grid
Anew = Anan<1;
gtest1 = g2D.xs{1}.*Anew; %removed not-allowed points from grid
gtest2 = g2D.xs{2}.*Anew;


gtest1_rem = g2D.xs{1}.*Anan; %removed not-allowed points from grid
gtest2_rem = g2D.xs{2}.*Anan;


% Plot feasible/infeasible states
figure(1)
ax1 = subplot(1,2,1);
suptitle('any velocity')


%get rid of all the 0's
ang1 = gtest1(gtest1~=0);%-pi/15;
ang2 = gtest2(gtest2~=0);%-pi/80;


%what the angles are relative to y axis
knee = ang1;
hip = ang2+ang1;

for j = 1:length(ang1); %for each state
    
angles=[knee(j) hip(j)]; %pull out the associated angles

x=zeros(1,4);
y=zeros(1,4);

L = [0 L0 L1 R2]; %length between each point
y(2)=L(2); %draw a straight line up from angle to knee

for i = 3:4 %find hip and head position
  x(i)=x(i-1) -  L(i)*sin(angles(i-2));
  y(i)=y(i-1) + L(i)*cos(angles(i-2));
end

plot(x,y,'Linewidth',3)
hold on
axis equal
xlabel('Feasible States')
end

ax2 = subplot(1,2,2);

%get rid of all the 0's
ang1 = gtest1_rem(gtest1_rem~=0);%-pi/15;
ang2 = gtest2_rem(gtest2_rem~=0);%-pi/80;

%what the angles are relative to y axis
knee = ang1;
hip = ang2+ang1;

for j = 1:length(ang1); %for each state
    
angles=[knee(j) hip(j)]; %pull out the associated angles

x=zeros(1,4);
y=zeros(1,4);

L = [0 L0 L1 R2]; %length between each point
y(2)=L(2); %draw a straight line up from angle to knee

for i = 3:4 %find hip and head position
  x(i)=x(i-1) -  L(i)*sin(angles(i-2));
  y(i)=y(i-1) + L(i)*cos(angles(i-2));
end

plot(x,y,'Linewidth',3)
hold on
end
xlabel('Infeasible States')
linkaxes([ax1,ax2],'xy')

elseif strcmp(view,'stand_conditions')
%%%NOTE!%%%
% if you want to see 2D through all possible velocities, use 'min'
% if you want to see 2D at a specific velocity, use [vKnee vHip], where 
% vKnee vHip are the velocities you want to use

[g2D, data2D]=proj(g,data,[0 1 0 1],'min');
gtest1=g2D.xs{1}.*(data2D<0);
gtest2=g2D.xs{2}.*(data2D<0);

figure(1)
title('Standing Conditions')

%get rid of all the 0's
ang1 = gtest1(gtest1~=0);%-pi/15;
ang2 = gtest2(gtest2~=0);%-pi/80;

%what the angles are relative to y axis
knee = ang1;
hip = ang2+ang1;

for j = 1:length(ang1); %for each state
    
angles=[knee(j) hip(j)]; %pull out the associated angles

x=zeros(1,4);
y=zeros(1,4);

L = [0 L0 L1 R2]; %length between each point
y(2)=L(2); %draw a straight line up from angle to knee

for i = 3:4 %find hip and head position
  x(i)=x(i-1) -  L(i)*sin(angles(i-2));
  y(i)=y(i-1) + L(i)*cos(angles(i-2));
end

plot(x,y,'Linewidth',3)
hold on
axis equal
xlabel('Feasible States')
end
end
