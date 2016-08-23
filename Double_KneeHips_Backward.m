
%% Grid
grid_min = [-1.8091-pi/15, -0.3289-pi/15, 0.0013-pi/15, -4.9904-pi/15];
grid_max = [0.0247+pi/15, 4.4269+pi/15, 2.4655+pi/15, 1.8952+pi/15];
%grid_min = [-1.8091-.5; -3.289-.5; .0013-.5; -4.99-.5]; % Lower corner of computation domain
%grid_max = [.0247+.5; 4.4269+.5; 2.4655+.5; 1.8952+.5];    % Upper corner of computation domain
N = 41*ones(4,1);         % Number of grid points per dimension
%pdDims = [1,3];               % 1st, 3rd dimension is periodic
g = createGrid(grid_min, grid_max, N);
% Use "g = createGrid(grid_min, grid_max, N);" if there are no periodic
% state space dimensions

%% target set
%data = shapeRectangleByCorners(g, [-pi/15;-pi/4; -pi/80;-pi/4],[pi/80;pi/4; pi/15;pi/4]);
 
%data = shapeRectangleByCorners(g, [-pi/15;-pi/15; -pi/80;-pi/15],[pi/80;pi/15; pi/15;pi/15]);
data = shapeRectangleByCorners(g, [-pi/15;-.1; -pi/80;-.1],[pi/80;.1; pi/15;.1]); 
% also try shapeRectangleByCorners, shapeSphere, etc.

%% time vector
t0 = 0;
tMax = 1;
dt = 0.025;
tau = t0:dt:tMax;
% If intermediate results are not needed, use tau = [t0 tMax];

%% problem parameters
T1Max = 1;
T2Max = 1;
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
%extraArgs.visualize = true;
extraArgs.save_filename = 'Double_KneeHips_Backward_41_highAc';
extraArgs.saveFrequency = 100;
[data, tau, extraOuts] = HJIPDE_solve( ...
  data, tau, schemeData, 'zero');

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