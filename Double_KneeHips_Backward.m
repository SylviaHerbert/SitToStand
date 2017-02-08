function [data, g, tau] = Double_KneeHips_Backward(gpoints, accuracy, tMax, T1Max, T1Min, T2Max, T2Min, TAMax, TAMin)
%% Grid
grid_min = [-0.29, -6.31, -2.67, -7.57];
grid_max = [1.89, 4.51, 0.15, 8.91];
N = gpoints*ones(4,1);         % Number of grid points per dimension, default to 41
g = createGrid(grid_min, grid_max, N);


%Want to drop # of test points from 8 to 6? use trim.
trim = 1;
%% target set
%data = shapeRectangleByCorners(g, [-pi/15;-pi/4; -pi/80;-pi/4],[pi/80;pi/4; pi/15;pi/4]);
 
%data = shapeRectangleByCorners(g, [-pi/15;-pi/15; -pi/80;-pi/15],[pi/80;pi/15; pi/15;pi/15]);

max_v = (pi/8)/2;       % allowing for some sway
standing_min = [-pi/15, -max_v, -pi/15, -max_v];
standing_max = [pi/15, max_v, 0.15, max_v];

data = shapeRectangleByCorners(g, standing_min, standing_max);
% also try shapeRectangleByCorners, shapeSphere, etc.

%% time vector
t0 = 0;

if nargin < 3
tMax = 2;
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

if nargin < 4
  T1Max = 107;
T1Min = -T1Max;
T2Max = 87;
T2Min = -60;
TAMax = 68;
TAMin = -50;
end

alpha = .1;
T1Max = alpha*T1Max;
T1Min = alpha*T1Min;
T2Max = alpha*T2Max;
T2Min = alpha*T2Min;
TAMax = alpha*TAMax;
TAMin = alpha*TAMin;
%% Pack problem parameters
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

% ----- System dynamics are specified here -----
schemeData.hamFunc = @pendulum4DHam;
schemeData.partialFunc = @pendulum4Dpartial;
%% solve
%extraArgs.visualize = true;
%extraArgs.save_filename = 'Double_KneeHips_Backward_41_highAc';
%extraArgs.saveFrequency = 100;
extraArgs.stopInit = [pi/2 0 -pi/2 0];
[data, tau] = HJIPDE_solve( ...
  data, tau, schemeData, 'zero', extraArgs);
%data = min(data,[],5);
end
% %% Ignore everything after this; using it to troubleshoot
% visualize = 1;
% schemeData.dissFunc = @artificialDissipationGLF;
% % Set up spatial approximation scheme.
% schemeFunc = @termLaxFriedrichs;
% accuracy = 'high';
% % Set up time approximation scheme.
% integratorOptions = odeCFLset('factorCFL', 0.5, 'stats', 'off');
% 
% switch(accuracy)
%   case 'low'
%     schemeData.derivFunc = @upwindFirstFirst;
%     integratorFunc = @odeCFL1;
%   case 'medium'
%     schemeData.derivFunc = @upwindFirstENO2;
%     integratorFunc = @odeCFL2;
%   case 'high'
%     schemeData.derivFunc = @upwindFirstENO3;
%     integratorFunc = @odeCFL3;
%   case 'veryHigh'
%     schemeData.derivFunc = @upwindFirstWENO5;
%     integratorFunc = @odeCFL3;
%   otherwise
%     error('Unknown accuracy level %s', accuracy);
% end
% integratorOptions = odeCFLset(integratorOptions, 'singleStep', 'on');
% 
% if visualize
%   figure(3);
%   clf
%   
%   [gKnee,dataKnee]=proj(g,data,[0 0 1 1]);
%   [~, h1] = contour(gKnee.xs{1}, gKnee.xs{2}, dataKnee, [0 0],'r');
%   hold on
%   xlabel('\theta Knee')
%   ylabel('\omega Knee')
%   drawnow;
%   
%   figure(4);
%   clf
%   [gHip,dataHip]=proj(g,data,[1 1 0 0]);
%   [~, h2] = contour(gHip.xs{1}, gHip.xs{2}, dataHip, [0 0],'r');
%   hold on
%   xlabel('\theta Hip')
%   ylabel('\omega Hip')
%   drawnow;
%   
% %   figure(3);
% %   clf
% %   [g3,data3]=proj(g,data,[0 1 0 1]);
% %   [~, h3] = contour(g3.xs{1}, g3.xs{2}, data3, [0 0],'r');
% %   hold on
% %   xlabel('\theta Knee')
% %   ylabel('\theta Hip')
% %   drawnow;
% %   
% %     figure(4);
% %   clf
% %   [g4,data4]=proj(g,data,[1 0 1 0]);
% %   [~, h4] = contour(g4.xs{1}, g4.xs{2}, data4, [0 0],'r');
% %   hold on
% %   xlabel('\omega Knee')
% %   ylabel('\omega Hip')
% %   drawnow;
% %   
% %       figure(5);
% %   clf
% %   [g5,data5]=proj(g,data,[1 0 0 1]);
% %   [~, h5] = contour(g5.xs{1}, g5.xs{2}, data5, [0 0],'r');
% %   hold on
% %   xlabel('\omega Knee')
% %   ylabel('\theta Hip')
% %   drawnow;
% %   
% %       figure(6);
% %   clf
% %   [g6,data6]=proj(g,data,[0 1 1 0]);
% %   [~, h6] = contour(g6.xs{1}, g6.xs{2}, data6, [0 0],'r');
% %   hold on
% %   xlabel('\theta Knee')
% %   ylabel('\omega Hip')
% %   drawnow;
% end
% 
% tNow = 0;
%  tau = tNow;
% 
% % How close (relative) do we need to get to tMax to be considered finished?
% small = 100 * eps;
% while(tMax - tNow > small * tMax)
%   % How far to step?
%   tSpan = [tNow, tMax];
%   
%     y0 = data(:,:,:,:,end);
%     y0 = y0(:);
%     [t, y] = feval(integratorFunc, schemeFunc, tSpan, y0,...
%       integratorOptions, schemeData);
%     %data = reshape(y,g.shape);
%         dataNew = reshape(y,g.shape);
%     data = min(data,dataNew);
% %     dataNew = min(data(:,:,:,:,end),dataNew);
% %     data = cat(5, data, dataNew);
%   
%     tNow = t(end);
%   tau = cat(1, tau, tNow);
%   
%   % Create new visualization.
%   if visualize
%     figure(3)
%     clf
%     [~, dataKnee]=proj(g,data,[0 0 1 1]);
%     %clf
%     [~, h1] = contour(gKnee.xs{1}, gKnee.xs{2}, dataKnee, [0 0],'r');
%     %h1.ZData = dataKnee;
%     title(['t=' num2str(tNow)])
%     drawnow;
%     
%     figure(4)
%     clf
%     [~, dataHip]=proj(g,data,[1 1 0 0]);
%     %clf
%     [~, h2] = contour(gHip.xs{1}, gHip.xs{2}, dataHip, [0 0],'r');
%     %h2.ZData = dataHips;
%     title(['t=' num2str(tNow)])
%     drawnow;
%     
% %   figure(3);
% %   [~,data3]=proj(g,data,[0 1 0 1]);
% %   [~, h3] = contour(g3.xs{1}, g3.xs{2}, data3, [0 0],'r');
% % title(['t=' num2str(tNow)])
% %   drawnow;
% %   
% %     figure(4);
% %   [~,data4]=proj(g,data,[1 0 1 0]);
% %   [~, h4] = contour(g4.xs{1}, g4.xs{2}, data4, [0 0],'r');
% % title(['t=' num2str(tNow)])
% %   drawnow;
% %   
% %       figure(5);
% %   [~,data5]=proj(g,data,[1 0 0 1]);
% %   [~, h5] = contour(g5.xs{1}, g5.xs{2}, data5, [0 0],'r');
% % title(['t=' num2str(tNow)])
% %   drawnow;
% %   
% %       figure(6);
% %   [~,data6]=proj(g,data,[0 1 1 0]);
% %   [~, h6] = contour(g6.xs{1}, g6.xs{2}, data6, [0 0],'r');
% % title(['t=' num2str(tNow)])
% %   drawnow;
%   end
% end
