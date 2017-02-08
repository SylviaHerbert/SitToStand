function simulationTest(schemeData, view)
%view options:
% allowed_states
% sitting_conditions
% standing_conditions
if nargin <2
  view = 'allowed_states';
end

% slice:
%'max' = union
%'min' = intersection
%[v1 v2] = slice at velocities
slice = 'max';

if nargin <1
  
  %% set up Grid, Initial Data
  
  %how many grid points we want in all directions
  gpoints = 45;
  N = gpoints*ones(4,1);         % Number of grid points per dimension, default to 41
  
  %bounds on all possible grid points
  grid_max = -[-1.8091-pi/15, -0.3289-pi/15, 0.0013-pi/15, -4.9904-pi/15];
  grid_min = -[0.0247+pi/15, 4.4269+pi/15, 2.4655+pi/15, 1.8952+pi/15];
  
  g = createGrid(grid_min, grid_max, N);
  
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
  trim = 1;
  schemeData = testAnkle(schemeData,trim);
end

%set what counts as "standing"
%data = shapeRectangleByCorners(g, lowerlimit, upperlimit)
max_v = pi/8;       % allowing for some sway
standing_min = [-pi/15, -max_v, -pi/15, -max_v]; 
standing_max = [pi/15, max_v, 0.15, max_v];
data = shapeRectangleByCorners(schemeData.grid, standing_min, standing_max);

%% finding qualified points within constraints
tau1Mins = schemeData.tau1test(:,:,:,:,1); %pull out grid of points that have feasible points to be tested
tau1Maxs = tau1Mins;
tau2Mins = schemeData.tau2test(:,:,:,:,1);
tau2Maxs = tau2Mins;

%b = 1;
%for each of the 6 possible points
for m = 2:length(schemeData.tau1test(1,1,1,1,:))
  tau1Mins = min(tau1Mins,schemeData.tau1test(:,:,:,:,m));
  tau1Maxs = max(tau1Maxs,schemeData.tau1test(:,:,:,:,m));
  tau2Mins = min(tau2Mins,schemeData.tau2test(:,:,:,:,m));
  tau2Maxs = max(tau2Maxs,schemeData.tau2test(:,:,:,:,m));
  %find all the states the have allowable torques
  %by taking the min across all of these points we're throwing out any
  %points that only have "NaN" as possible torques
end
tau1Minstemp = (tau1Mins==0);
tau1Maxstemp = (tau1Maxs==0);
tau1NoControl = ((tau1Minstemp+tau1Maxstemp) == 2);

tau2Minstemp = (tau2Mins==0);
tau2Maxstemp = (tau2Maxs==0);
tau2NoControl = ((tau2Minstemp+tau2Maxstemp) == 2);

NoControl = ((tau1NoControl+tau2NoControl)==2); %points where min and max for everything is 0
Control = 1-NoControl;

%% Project data through to 2D so we can view it

if strcmp(view,'allowed_states')
  %%%NOTE!%%%
  % if you want to see 2D through all possible velocities, use 'max'
  % if you want to see 2D through intersection of velocities, use 'min'
  % if you want to see 2D at a specific velocity, use [vKnee vHip], where
  % vKnee vHip are the velocities you want to use
  
  [g2D, Control2D]=proj(schemeData.grid,Control,[0 1 0 1],slice); %project through to 2D
  
  % we care only about the points on this 2D grid that have feasible torques.
  % we're going to do some awkward manipulation to make sure we only have
  % those:
  gtest1 = g2D.xs{1}.*Control2D; %removed not-allowed points from grid
  gtest2 = g2D.xs{2}.*Control2D;
  
  
  gtest1_rem = g2D.xs{1}.*(1-Control2D); %removed not-allowed points from grid
  gtest2_rem = g2D.xs{2}.*(1-Control2D);
  
  
  % Plot feasible/infeasible states
  figure(1)
  clf
  ax1 = subplot(1,2,1);
  if strcmp(slice,'min')
    suptitle('every velocity (intersection)')
  elseif strcmp(slice,'max')
    suptitle('any velocity (union)')
  else
    suptitle(['vKnee = ' num2str(slice(1)) ', vHip = ' num2str(slice(2))])
  end
  
  
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
    
    L = [0 schemeData.L0 schemeData.L1 schemeData.R2]; %length between each point
    y(2)=L(2); %draw a straight line up from angle to knee
    
    for i = 3:4 %find hip and head position
      x(i)=x(i-1) -  L(i)*sin(angles(i-2));
      y(i)=y(i-1) + L(i)*cos(angles(i-2));
    end
    
    plot(x,y,'Linewidth',3)
    hold on
    axis([-1 1 0 1.5])
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
    
    L = [0 schemeData.L0 schemeData.L1 schemeData.R2]; %length between each point
    y(2)=L(2); %draw a straight line up from angle to knee
    
    for i = 3:4 %find hip and head position
      x(i)=x(i-1) -  L(i)*sin(angles(i-2));
      y(i)=y(i-1) + L(i)*cos(angles(i-2));
    end
    
    plot(x,y,'Linewidth',3)
    hold on
  end
  xlabel('Infeasible States')
  axis([-1 1 0 1.5])
  %linkaxes([ax1,ax2],'xy')
  
elseif strcmp(view,'sitting_conditions')
  %%%NOTE!%%%
  % if you want to see 2D through all possible velocities, use 'max'
  % if you want to see 2D through intersection of velocities, use 'min'
  % if you want to see 2D at a specific velocity, use [vKnee vHip], where
  % vKnee vHip are the velocities you want to use
  
  [g2D, Control2D]=proj(schemeData.grid,Control,[0 1 0 1],slice); %project through to 2D
  [g1D, Control1D]=proj(g2D,Control2D,[1 0],pi/2);
  
  % we care only about the points on this 2D grid that have feasible torques.
  % we're going to do some awkward manipulation to make sure we only have
  % those:
  gtest = g1D.xs{1}.*Control1D; %removed not-allowed points from grid
  
  
  gtest_rem = g1D.xs{1}.*(1-Control1D); %removed not-allowed points from grid
  
  
  % Plot feasible/infeasible states
  figure(1)
  clf
  ax1 = subplot(1,2,1);
  if strcmp(slice,'min')
    suptitle('sitting at every velocity (intersection)')
  elseif strcmp(slice,'max')
    suptitle('sitting at any velocity (union)')
  else
    suptitle(['sitting at vKnee = ' num2str(slice(1)) ', vHip = ' num2str(slice(2))])
  end
  
  
  %get rid of all the 0's
  ang2 = gtest(gtest~=0);
  ang1 = ones(length(ang2),1).*(pi/2);
  
  
  %what the angles are relative to y axis
  knee = ang1;
  hip = ang2+ang1;
  
  for j = 1:length(ang1); %for each state
    
    angles=[knee(j) hip(j)]; %pull out the associated angles
    
    x=zeros(1,4);
    y=zeros(1,4);
    
    L = [0 schemeData.L0 schemeData.L1 schemeData.R2]; %length between each point
    y(2)=L(2); %draw a straight line up from angle to knee
    
    for i = 3:4 %find hip and head position
      x(i)=x(i-1) -  L(i)*sin(angles(i-2));
      y(i)=y(i-1) + L(i)*cos(angles(i-2));
    end
    
    plot(x,y,'Linewidth',3)
    hold on
    axis([-1 1 0 1.5])
    %axis equal
    xlabel('Feasible States')
  end
  
  ax2 = subplot(1,2,2);
  
  %get rid of all the 0's
  ang2 = gtest_rem(gtest_rem~=0);%-pi/15;
  ang1 = ones(length(ang2),1).*(pi/2);
  
  %what the angles are relative to y axis
  knee = ang1;
  hip = ang2+ang1;
  
  for j = 1:length(ang1); %for each state
    
    angles=[knee(j) hip(j)]; %pull out the associated angles
    
    x=zeros(1,4);
    y=zeros(1,4);
    
    L = [0 schemeData.L0 schemeData.L1 schemeData.R2]; %length between each point
    y(2)=L(2); %draw a straight line up from angle to knee
    
    for i = 3:4 %find hip and head position
      x(i)=x(i-1) -  L(i)*sin(angles(i-2));
      y(i)=y(i-1) + L(i)*cos(angles(i-2));
    end
    
    plot(x,y,'Linewidth',3)
    hold on
  end
  xlabel('Infeasible States')
  axis([-1 1 0 1.5])
  %linkaxes([ax1,ax2],'xy')
  
elseif strcmp(view,'standing_conditions')
  %%%NOTE!%%%
  % if you want to see 2D through all possible velocities, use 'min'
  % if you want to see 2D at a specific velocity, use [vKnee vHip], where
  % vKnee vHip are the velocities you want to use
  
  [g2D, data2D]=proj(schemeData.grid,data,[0 1 0 1],'min');
  gtest1=g2D.xs{1}.*(data2D<=0);
  gtest2=g2D.xs{2}.*(data2D<=0);
  
  gtest1_rem = g2D.xs{1}.*(data2D>0); %removed not-allowed points from grid
  gtest2_rem = g2D.xs{2}.*(data2D>0);
  
  figure(1)
  clf
  ax2 = subplot(1,2,1);
  axis([-1 1 0 1.5])
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
    
    L = [0 schemeData.L0 schemeData.L1 schemeData.R2]; %length between each point
    y(2)=L(2); %draw a straight line up from angle to knee
    
    for i = 3:4 %find hip and head position
      x(i)=x(i-1) -  L(i)*sin(angles(i-2));
      y(i)=y(i-1) + L(i)*cos(angles(i-2));
    end
    
    plot(x,y,'Linewidth',3)
    hold on
    axis([-1 1 0 1.5])
    xlabel('Feasible States')
  end
  
  
    ax2 = subplot(1,2,2);
  
  %get rid of all the 0's
  ang1 = gtest1_rem(gtest1_rem~=0);%-pi/15;
  ang2 = gtest2_rem(gtest2_rem~=0);%-pi/15;
  
  %what the angles are relative to y axis
  knee = ang1;
  hip = ang2+ang1;
  
  for j = 1:length(ang1); %for each state
    
    angles=[knee(j) hip(j)]; %pull out the associated angles
    
    x=zeros(1,4);
    y=zeros(1,4);
    
    L = [0 schemeData.L0 schemeData.L1 schemeData.R2]; %length between each point
    y(2)=L(2); %draw a straight line up from angle to knee
    
    for i = 3:4 %find hip and head position
      x(i)=x(i-1) -  L(i)*sin(angles(i-2));
      y(i)=y(i-1) + L(i)*cos(angles(i-2));
    end
    
    plot(x,y,'Linewidth',3)
    hold on
  end
  xlabel('Infeasible States')
  %axis equal
  axis([-1 1 0 1.5])
  %linkaxes([ax1,ax2],'xy')
end
end