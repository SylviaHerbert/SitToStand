function body=simulateTrajectories(g,data,tau,z0,schemeData)
uMode = 'min';
dt = .01;


if nargin < 5
  height = 1.72;
mass = 62;
schemeData.M1 = 2*(0.1416*mass);       % mass of thighs 
schemeData.M2 = (.0694 +.4346)*mass;    % mass of head-arms-trunk
schemeData.L0 = .25*height;          % length of segment (shank)
schemeData.L1 = .26*height;          % length of segment .4(thigh)
schemeData.R1 = .43*schemeData.L1;          % position of COM along segment (thigh)
schemeData.R2 = .6*.4*height;          % position of COM along segment (head-arms-trunk)
schemeData.L2 = .3*height;
schemeData.grav = 9.81;

T1Max = 107;
T1Min = -T1Max;
T2Max = 87;
T2Min = -60;
TAMax = 68;
TAMin = -50;


alpha = 1;
schemeData.T1Max = alpha*T1Max;
schemeData.T1Min = alpha*T1Min;
schemeData.T2Max = alpha*T2Max;
schemeData.T2Min = alpha*T2Min;
schemeData.TAMax = alpha*TAMax;
schemeData.TAMin = alpha*TAMin;
end

R1 = schemeData.R1;
R2 = schemeData.R2;
M1 = schemeData.M1;
M2 = schemeData.M2;
L1 = schemeData.L1;
L0 = schemeData.L0;
grav = 9.81;
T1Max = schemeData.T1Max;
T1Min = schemeData.T1Min;
T2Max = schemeData.T2Max;
T2Min = schemeData.T2Min;
TAMax = schemeData.TAMax;
TAMin = schemeData.TAMin;

dims = 1:4;

body = STS4D(z0, R1, R2, M1, M2, L1, L0, grav,...
        T1Max, T1Min, T2Max, T2Min, TAMax, TAMin, dims);
%% Initial figure
f=figure(1);
clf
[x,y]=transferCoordinates(z0,body);
h=plot(x,y,'Linewidth',3);
axis([-1 1 0 1.5])
%% Find all Gradients
% for i = 1:length(tau)
%   deriv{i} = computeGradients(g,data(:,:,:,:,i));
% end

%% Loop over every time stamp

for i = 1:length(data(1,1,1,1,:))
  deriv{i} = computeGradients(g,data(:,:,:,:,i));
end

for i = 2:length(tau)
% Find gradient at this state
p = eval_u(g,deriv{i},body.x);

% Find Optimal Control
% at this state, what are the allowed controls?
% Find the best of those controls

%[body] = tauPoints(body);

uOpt = body.optCtrl(dt,body.x,p,uMode);

% Update State
body.updateState(uOpt, dt);
z = body.x;

% Plot
[x,y]=transferCoordinates(z,body);
h.XData = x;
h.YData = y;
title(num2str(tau(i)));
pause

max_v = (pi/8);       % allowing for some sway
standing_min = [-pi/15, -max_v, -pi/15, -max_v]';
standing_max = [pi/15, max_v, 0.15, max_v]';
upper = (z <= standing_max);
lower = (z >= standing_min);

if sum(upper)==4 && sum(lower)==4
  break
end

end
end

function [obj] = tauPoints(obj)

x = obj.x;

%g=schemeData.grid; % Grid MUST be specified!
tau1Bound = [obj.T1Min, obj.T1Max];
tau2Bound = [obj.T2Min, obj.T2Max];
tauAnkleMin = obj.TAMin;
tauAnkleMax = obj.TAMax;
R1=obj.R1;
R2=obj.R2;
M1=obj.M1;
M2=obj.M2;
L1=obj.L1;
L0=obj.L0;
grav = 9.81;
%dx(1) = x(2);

%dx(2) = tau1.*(tau1num1/denom1) + tau2.*(tau2num1/denom1) + (num1/denom1);
denom1 = (L1.^2.*M2.*R2 + M1.*R1.^2.*R2 - ...
  L1.^2.*M2.*R2.*cos(x(3)).^2);
num1 = grav.*(M1.*R1.*R2 + M2.*L1.*R2).*sin(x(1)) + ...
  (x(2) + x(4)).^2.*L1.*M2.*R2.^2.*sin(x(3)) - ...
  grav.*M2.*L1.*R2.*sin(x(1)+x(3)).*cos(x(3)) + ...
  M2.*L1.^2.*R2.*x(2).^2.*cos(x(3)).*sin(x(3));
tau1num1 = R2;
tau2num1 = -(R2+L1.*cos(x(3)));

%dx(3) = x(4);

%dx(4) = tau1.*(tau1num2/denom2) + tau2.*(tau2num2/denom2)  + (num2/denom2);
num2 = -((M2.^2.*R2.^2.*L1.*grav + M1.*M2.*R1.*R2.^2.*grav).*sin(x(1)) + ...
  (-M2.^2.*R2.*L1.^2.*grav - M1.*M2.*R1.^2.*R2.*grav).*sin(x(1) + x(3)) + ...
  ((M2.^2.*R2.*L1.^3+M1.*M2.*R1.^2.*R2.*L1).*x(2).^2+ ...
  (M2.^2.*R2.^3.*L1).*(x(2)+x(4)).^2).*sin(x(3)) + ...
  (M2.^2.*R2.*L1.^2.*grav + M1.*M2.*R1.*R2.*L1.*grav).*cos(x(3)).*sin(x(1)) + ...
  (M2.^2.*R2.^2.*L1.^2.*(2.*x(2).^2 + 2.*x(2).*x(4) + x(4).^2)).*cos(x(3)).*sin(x(3)) - ...
  M2.^2.*R2.^2.*L1.*grav.*sin(x(1) + x(3)).*cos(x(3)));
denom2 = (M2.*R2.^2.*(M1.*R1.^2 + M2.*L1.^2 - M2.*L1.^2.*cos(x(3)).^2));
tau1num2 = -(M2.*R2.^2 + M2.*R2.*L1.*cos(x(3)));
tau2num2= -(-M1.*R1.^2 - M2.*R2.^2 - M2.*L1.^2 - 2.*M2.*R2.*L1.*cos(x(3)));


%%

%ankle torque constraint:  tau_min < tao0 < tau_max


tau0num1 = (L1.^2.*M2 + M1.*R1.^2 + M2.*R2.^2 + L0.*M2.*R2.*cos(x(1) ...
  + x(3)) + L0.*L1.*M2.*cos(x(1)) + L0.*M1.*R1.*cos(x(1)) + ...
  L1.*M2.*R2.*cos(x(3)));
tau0num2 = (M2.*R2.^2 + L0.*M2.*R2.*cos(x(1) + x(3)) + L1.*M2.*R2.*cos(x(3)));
tau0num3 = - M2.*R2.*grav.*sin(x(1) + x(3))...
  - (L1.*M2.*grav + M1.*R1.*grav).*sin(x(1)) + ...
  - L0.*M2.*R2.*x(2).^2.*sin(x(1) + x(3)) - L0.*L1.*M2.*x(2).^2.*sin(x(1)) ...
  - L0.*M2.*R2.*x(4).^2.*sin(x(1) + x(3)) ...
  - L0.*M1.*R1.*x(2).^2.*sin(x(1)) - L1.*M2.*R2.*x(4).^2.*sin(x(3)) - ...
  L0.*M2.*R2.*0.*x(2).*sin(x(1) + x(3)) - L0.*M2.*R2.*0.*x(4).*sin(x(1) ...
  + x(3)) - L0.*M2.*R2.*x(2).*x(4).*sin(x(1) + x(3)) - ...
  (L0.*L1.*M2.*0 + L0.*M1.*R1.*0).*x(2).*sin(x(1)) ...
  - L1.*M2.*R2.*0.*x(4).*sin(x(3)) - L1.*M2.*R2.*x(2).*x(4).*sin(x(3));

tau1multiplier = tau0num1.*(tau1num1./denom1)+tau0num2.*(tau1num2./denom2);
tau2multiplier = tau0num1.*(tau2num1./denom1)+tau0num2.*(tau2num2./denom2);
extraTerms = (num1./denom1).*tau0num1 + (num2./denom2).*tau0num2 + tau0num3;

tau1Test = [];
tau2Test = [];
for i = 1:2
  tau1 = tau1Bound(i);
  
  tau2Lim1 = (tauAnkleMin - tau1.*tau1multiplier - extraTerms)./tau2multiplier;
  tau2Lim2 = (tauAnkleMax - tau1.*tau1multiplier - extraTerms)./tau2multiplier;
  
  tau2MinTemp = min(tau2Lim1, tau2Lim2);
  tau2MaxTemp = max(tau2Lim1, tau2Lim2);
  
  tau2Min = max(tau2MinTemp, tau2Bound(1));
  tau2Max = min(tau2MaxTemp, tau2Bound(2));
  
  if tau2Max >= tau2Min
    tau1Test = [tau1Test, tau1, tau1];
    tau2Test = [tau2Test, tau2Min, tau2Max];
  end
end
for i = 1:2
  tau2 = tau2Bound(i);
  
  tau1Lim1 = (tauAnkleMin - tau2.*tau2multiplier - extraTerms)./tau1multiplier;
  tau1Lim2 = (tauAnkleMax - tau2.*tau2multiplier - extraTerms)./tau1multiplier;
  
  tau1MinTemp = min(tau1Lim1, tau1Lim2);
  tau1MaxTemp = max(tau1Lim1, tau1Lim2);
  
  tau1Min = max(tau1MinTemp, tau1Bound(1));
  tau1Max = min(tau1MaxTemp, tau1Bound(2));
  
  if tau1Max >= tau1Min
    tau1Test = [tau1Test, tau1Min, tau1Max];
    tau2Test = [tau2Test, tau2, tau2];
  end
  obj.tau1Test = tau1Test;
  obj.tau2Test = tau2Test;
end
end

function [x,y]=transferCoordinates(z,body)
x=zeros(1,4);
y=zeros(1,4);
L = [0 body.L0 body.L1 body.R2]; %length between each point
y(2)=L(2); %draw a straight line up from angle to knee
x(3) = x(2)-L(3)*sin(z(1));
y(3) = y(2) + L(3)*cos(z(1));
x(4) = x(3)-L(4)*sin(z(1)+z(3));
y(4) = y(3)+L(4)*cos(z(1)+z(3));
end