%get the grid and the initial data 

gpoints=45;
grid_min = [-1.8091-pi/15, -0.3289-pi/15, 0.0013-pi/15, -4.9904-pi/15];
grid_max = [0.0247+pi/15, 4.4269+pi/15, 2.4655+pi/15, 1.8952+pi/15];
%grid_min = [-1.8091-.5; -3.289-.5; .0013-.5; -4.99-.5]; % Lower corner of computation domain
%grid_max = [.0247+.5; 4.4269+.5; 2.4655+.5; 1.8952+.5];    % Upper corner of computation domain
N = gpoints*ones(4,1);         % Number of grid points per dimension, default to 41
%pdDims = [1,3];               % 1st, 3rd dimension is periodic
g = createGrid(grid_min, grid_max, N);

data = shapeRectangleByCorners(g, [-pi/15;-.1; -pi/80;-.1],[pi/80;.1; pi/15;.1]); 

height = 1.72;
mass = 62;
M1 = 2*(0.1416*mass);       % mass of thighs 
M2 = (.0694 +.4346)*mass;    % mass of head-arms-trunk
L0 = .25*height;          % length of segment (shank)
L1 = .26*height;          % length of segment .4(thigh)
R1 = .43*L1;          % position of COM along segment (thigh)
R2 = .6*.4*height;  
  T1Max = 107;
  T1Min = -T1Max;
  T2Max = 87;
  T2Min = -60;
  TAMax = 68;
  TAMin = -50;
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
schemeData = testAnkle(schemeData);


%% finding qualified points within constraints
A = schemeData.tau1test{1};
for m = 2:8
  A = min(A,schemeData.tau1test{m}); %find all the states the have allowable torques
end
[g2D, A2D]=proj(g,A,[0 1 0 1],'min'); %project through to 2D
Anan = isnan(A2D); %find points that are not allowed on 2D grid
Anew = Anan<1;
gtest1 = g2D.xs{1}.*Anew; %removed not-allowed points form grid
gtest2 = g2D.xs{2}.*Anew;
%% test

% if solving for standing positions, un-comment this code
%[g2D, data2D]=proj(g,data,[0 1 0 1],'min');
%gtest1=g2D.xs{1}.*(data2D<0);
%gtest2=g2D.xs{2}.*(data2D<0);

%put all the angles together
ang1 = gtest1(gtest1~=0);%-pi/15;
ang2 = gtest2(gtest2~=0);%-pi/80;

%what the angles are relative to y axis
knee = ang1;
hip = ang2+ang1;

clf
for j = 1:length(ang1); %for each allowed state
angles=[knee(j) hip(j)];

height = 1.72;
mass = 62;
M1 = 2*(0.1416*mass);       % mass of thighs 
M2 = (.0694 +.4346)*mass;    % mass of head-arms-trunk
L0 = .25*height;          % length of segment (shank)
L1 = .26*height;          % length of segment .4(thigh)
L2 = .49*height; %(torso)
R1 = .43*L1;          % position of COM along segment (thigh)
R2 = .6*.4*height;

x=zeros(1,4);
y=zeros(1,4);

L = [0 L0 L1 L2];
y(2)=L(2);

for i = 3:4 %find hip and head position
  x(i)=x(i-1)+ L(i)*sin(angles(i-2));
  y(i)=y(i-1)+ L(i)*cos(angles(i-2));
end

plot(x,y,'Linewidth',3)
axis([-height height 0 1.5*height])
hold on
end