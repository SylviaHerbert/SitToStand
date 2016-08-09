
%% Grid
grid_min = [-1.8091-.5; -3.289-.5; .0013-.5; -4.99-.5]; % Lower corner of computation domain
grid_max = [.0247+.5; 4.4269+.5; 2.4655+.5; 1.8952+.5];    % Upper corner of computation domain
N = [61; 61; 61; 61];         % Number of grid points per dimension
%pdDims = [1,3];               % 1st, 3rd dimension is periodic
g = createGrid(grid_min, grid_max, N);
% Use "g = createGrid(grid_min, grid_max, N);" if there are no periodic
% state space dimensions

%% target set

data = shapeRectangleByCorners(g, [-pi/15;-.1; -pi/80;-.1],[pi/80;.1; pi/15;.1]);

% also try shapeRectangleByCorners, shapeSphere, etc.

%% time vector
t0 = 0;
tMax = 1;
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
schemeData.accuracy = 'medium';

% ----- System dynamics are specified here -----
schemeData.hamFunc = @pendulum4DHam;
schemeData.partialFunc = @pendulum4Dpartial;
%% solve
extraArgs.visualize = true;
extraArgs.save_filename = 'Double_KneeHips_Backward_61';
extraArgs.saveFrequency = 100;
[data, tau, extraOuts] = HJIPDE_solve( ...
  data, tau, schemeData, 'zero');

% visualize = 1;
% schemeData.dissFunc = @artificialDissipationGLF;
% % Set up spatial approximation scheme.
% schemeFunc = @termLaxFriedrichs;
% accuracy = 'medium';
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
%   figure(5);
%   clf
%   
%   [gKnee,dataKnee]=proj(g,data,[0 0 1 1]);
%   [~, h1] = contour(gKnee.xs{1}, gKnee.xs{2}, dataKnee, [0 0],'r');
%   hold on
%   xlabel('\theta Knee')
%   ylabel('\omega Knee')
%   drawnow;
%   
%   figure(6);
%   clf
%   [gHip,dataHip]=proj(g,data,[1 1 0 0]);
%   [~, h2] = contour(gHip.xs{1}, gHip.xs{2}, dataHip, [0 0],'r');
%   hold on
%   xlabel('\theta Hip')
%   ylabel('\omega Hip')
%   drawnow;
% end
% 
% tNow = 0;
% tau = tNow;
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
%     dataNew = reshape(y,g.shape);
%     data = min(data,dataNew);
% %     dataNew = min(data(:,:,:,:,end),dataNew);
% %     data = cat(5, data, dataNew);
%   
%     tNow = t(end);
%   tau = cat(1, tau, tNow);
%   
%   % Create new visualization.
%   if visualize
%     figure(5)
%     [~, dataKnee]=proj(g,data,[0 0 1 1]);
%     clf
%     [~, h1] = contour(gKnee.xs{1}, gKnee.xs{2}, dataKnee, [0 0],'r');
%     %h1.ZData = dataKnee;
%     title(['t=' num2str(tNow)])
%     drawnow;
%     
%     figure(6)
%     [~, dataHip]=proj(g,data,[1 1 0 0]);
%     clf
%     [~, h2] = contour(gHip.xs{1}, gHip.xs{2}, dataHip, [0 0],'r');
%     %h2.ZData = dataHips;
%     title(['t=' num2str(tNow)])
%     drawnow;
% 
%   end
% end