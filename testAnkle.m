function [schemeData] = testAnkle(schemeData, trim)
if nargin <1
grid_min = [-0.29, -6.31, -2.67, -7.57];
grid_max = [1.89, 4.51, 0.15, 8.91];
N = 20*ones(4,1);         % Number of grid points per dimension, default to 41
g = createGrid(grid_min, grid_max, N);
height = 1.72;
mass = 62;
M1 = 2*(0.1416*mass);       % mass of thighs 
M2 = (.0694 +.4346)*mass;    % mass of head-arms-trunk
L0 = .25*height;          % length of segment (shank)
L1 = .26*height;          % length of segment .4(thigh)
R1 = .43*L1;          % position of COM along segment (thigh)
R2 = .6*.4*height;          % position of COM along segment (head-arms-trunk)
L2 = .3*height;

  T1Max = 107;
T1Min = -T1Max;
T2Max = 87;
T2Min = -60;
TAMax = 68;
TAMin = -50;

alpha = .1;
T1Max = alpha*T1Max;
T1Min = alpha*T1Min;
T2Max = alpha*T2Max;
T2Min = alpha*T2Min;
TAMax = alpha*TAMax;
TAMin = alpha*TAMin;

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
schemeData.L2 = L2;
end

if nargin <2
  trim = 1;
end
%%% double pendulum - ankle torque constraint
%% parameters
g=schemeData.grid; % Grid MUST be specified!
tau1Bound = [schemeData.T1Min, schemeData.T1Max];
tau2Bound = [schemeData.T2Min, schemeData.T2Max];
tauAnkleMin = schemeData.TAMin;
tauAnkleMax = schemeData.TAMax;
R1=schemeData.R1;
R2=schemeData.R2;
M1=schemeData.M1;
M2=schemeData.M2;
L1=schemeData.L1;
L0=schemeData.L0;
grav = 9.81;

%% Setting up the Ankle Constraint Equation
x1 = g.xs{1};
x2 = g.xs{2};
x3 = g.xs{3};
x4 = g.xs{4};

%dx1 = x2;

%dx2 = tau1.*(tau1num1/denom1) + tau2.*(tau2num1/denom1) + (num1/denom1);
denom1 = (L1.^2.*M2.*R2 + M1.*R1.^2.*R2 - ...
  L1.^2.*M2.*R2.*cos(x3).^2);
num1 = grav.*(M1.*R1.*R2 + M2.*L1.*R2).*sin(x1) + ...
  (x2 + x4).^2.*L1.*M2.*R2.^2.*sin(x3) - ...
  grav.*M2.*L1.*R2.*sin(x1+x3).*cos(x3) + ...
  M2.*L1.^2.*R2.*x2.^2.*cos(x3).*sin(x3);
tau1num1 = R2;
tau2num1 = -(R2+L1.*cos(x3));

%dx3 = x4;

%dx4 = tau1.*(tau1num2/denom2) + tau2.*(tau2num2/denom2)  + (num2/denom2);
num2 = -((M2.^2.*R2.^2.*L1.*grav + M1.*M2.*R1.*R2.^2.*grav).*sin(x1) + ...
  (-M2.^2.*R2.*L1.^2.*grav - M1.*M2.*R1.^2.*R2.*grav).*sin(x1 + x3) + ...
  ((M2.^2.*R2.*L1.^3+M1.*M2.*R1.^2.*R2.*L1).*x2.^2+ ...
  (M2.^2.*R2.^3.*L1).*(x2+x4).^2).*sin(x3) + ...
  (M2.^2.*R2.*L1.^2.*grav + M1.*M2.*R1.*R2.*L1.*grav).*cos(x3).*sin(x1) + ...
  (M2.^2.*R2.^2.*L1.^2.*(2.*x2.^2 + 2.*x2.*x4 + x4.^2)).*cos(x3).*sin(x3) - ...
  M2.^2.*R2.^2.*L1.*grav.*sin(x1 + x3).*cos(x3));
denom2 = (M2.*R2.^2.*(M1.*R1.^2 + M2.*L1.^2 - M2.*L1.^2.*cos(x3).^2));
tau1num2 = -(M2.*R2.^2 + M2.*R2.*L1.*cos(x3));
tau2num2= -(-M1.*R1.^2 - M2.*R2.^2 - M2.*L1.^2 - 2.*M2.*R2.*L1.*cos(x3));


%%

%ankle torque constraint:  tau_min < tao0 < tau_max


tau0num1 = (L1.^2.*M2 + M1.*R1.^2 + M2.*R2.^2 + L0.*M2.*R2.*cos(x1 ...
  + x3) + L0.*L1.*M2.*cos(x1) + L0.*M1.*R1.*cos(x1) + ...
  L1.*M2.*R2.*cos(x3));
tau0num2 = (M2.*R2.^2 + L0.*M2.*R2.*cos(x1 + x3) + L1.*M2.*R2.*cos(x3));
tau0num3 = - M2.*R2.*grav.*sin(x1 + x3)...
  - (L1.*M2.*grav + M1.*R1.*grav).*sin(x1) + ...
  - L0.*M2.*R2.*x2.^2.*sin(x1 + x3) - L0.*L1.*M2.*x2.^2.*sin(x1) ...
  - L0.*M2.*R2.*x4.^2.*sin(x1 + x3) ...
  - L0.*M1.*R1.*x2.^2.*sin(x1) - L1.*M2.*R2.*x4.^2.*sin(x3) - ...
  L0.*M2.*R2.*0.*x2.*sin(x1 + x3) - L0.*M2.*R2.*0.*x4.*sin(x1 ...
  + x3) - L0.*M2.*R2.*x2.*x4.*sin(x1 + x3) - ...
  (L0.*L1.*M2.*0 + L0.*M1.*R1.*0).*x2.*sin(x1) ...
  - L1.*M2.*R2.*0.*x4.*sin(x3) - L1.*M2.*R2.*x2.*x4.*sin(x3);

%MATH!
% From Sarah's code:
% tau0 = dx2.*tau0num1 + dx4.*tau0num2 + tau0num3;

% plug in dx2, dx4 from above:
% tau0 = (tau1.*(tau1num1/denom1) - tau2.*(tau2num1/denom1) +
% (num1/denom1)).*tau0num1 +...
% (tau1.*(tau1num2/denom2) + tau2.*(tau2num2/denom2)  +
% (num2/denom2)).*tau0num2 + ...
% tau0num3;

% rearrange to separate tau1 and tau2:
% tau0 = tau1.*tau0num1.*(tau1num1/denom1) +tau1.*(tau1num2/denom2).*tau0num2 
% + tau2.*tau0num1.*(tau2num1/denom1) + tau2.*(tau2num2/denom2).*tau0num2 
% +(num1/denom1).*tau0num1 +(num2/denom2).*tau0num2 +tau0num3;


% simplify:
% tau0 = tau1.*(tau0num1.*(tau1num1./denom1)+(tau1num2./denom2).*tau0num2)...
% + tau2.*(tau0num1.*(tau2num1./denom1) + (tau2num2./denom2).*tau0num2) 
% +(num1./denom1).*tau0num1 +(num2./denom2).*tau0num2 +tau0num3;

tau1multiplier = tau0num1.*(tau1num1./denom1)+tau0num2.*(tau1num2./denom2);
tau2multiplier = tau0num1.*(tau2num1./denom1)+tau0num2.*(tau2num2./denom2);
extraTerms = (num1./denom1).*tau0num1 + (num2./denom2).*tau0num2 + tau0num3;

% tau0 = tau1.*(tau1multiplier) + tau2.*(tau2multiplier) + extraTerms;




%% Solving for test points

testPoints = {};
tau1TestPoints = [];
tau2TestPoints = [];
% First solve for tau1 = T1Min then tau1 = T1Max
for i = 1:2
  tau1 = tau1Bound(i);
  
  % Ankle Constraint that we're solving for:
  % tau0 = tau1.*(tau1multiplier) + tau2.*(tau2multiplier) + extraTerms;
  % Write this in terms of tau 2:
  %tau2 = (tau0 - tau1.*tau1multiplier - extraTerms)./tau2multiplier

  % set tau0 to be tauAnkleMin, then tauAnkleMax
  tau2Lim1 = (tauAnkleMin - tau1.*tau1multiplier - extraTerms)./tau2multiplier;
  tau2Lim2 = (tauAnkleMax - tau1.*tau1multiplier - extraTerms)./tau2multiplier;
 
  tau2MinTemp = min(tau2Lim1, tau2Lim2);
  tau2MaxTemp = max(tau2Lim1, tau2Lim2);
  
  %this gives us the min and max values of tau2 that the ankle constraint
  %will allow us.  Let's compare it to the nominal min and max values of
  %tau2 that we set.  We'll want to take the intersection between these two
  %ranges to get the range that both constraints allow.
  
  tau2Min = max(tau2MinTemp, tau2Bound(1));
  tau2Max = min(tau2MaxTemp, tau2Bound(2));
  
  %tau2Max must be greater than tau2Min for this to be a nonempty set
  A = tau2Max>=tau2Min; %cases when this is a nonempty set
  tau1 = tau1.*A;
  
  tau1TestPoints = cat(5,tau1TestPoints, tau1,tau1);
  tau2TestPoints = cat(5,tau2TestPoints, tau2Min.*A, tau2Max.*A);
  testPoints{end+1} = cat(5,tau1,tau2Min.*A);
  testPoints{end+1} = cat(5,tau1,tau2Max.*A);
end

%Now solve for tau2 = T2Min then tau2 = T2Max
for i = 1:2
  tau2 = tau2Bound(i);
  
  % Ankle Constraint that we're solving for:
  % tau0 = tau1.*(tau1multiplier) + tau2.*(tau2multiplier) + extraTerms;
  % Write this in terms of tau 1:
  % tau1 = (tau0 - tau2.*tau2multiplier - extraTerms)./tau1multiplier

  % set tau0 to be tauAnkleMin, then tauAnkleMax
  tau1Lim1 = (tauAnkleMin - tau2.*tau2multiplier - extraTerms)./tau1multiplier;
  tau1Lim2 = (tauAnkleMax - tau2.*tau2multiplier - extraTerms)./tau1multiplier;
 
  tau1MinTemp = min(tau1Lim1, tau1Lim2);
  tau1MaxTemp = max(tau1Lim1, tau1Lim2);
  
  %this gives us the min and max values of tau2 that the ankle constraint
  %will allow us.  Let's compare it to the nominal min and max values of
  %tau2 that we set.  We'll want to take the intersection between these two
  %ranges to get the range that both constraints allow.
  
  tau1Min = max(tau1MinTemp, tau1Bound(1));
  tau1Max = min(tau1MaxTemp, tau1Bound(2));
  
  %tau2Max must be greater than tau2Min for this to be a nonempty set
  A = tau1Max>=tau1Min; %cases when this is a nonempty set
  tau2 = tau2.*A;
  
  tau1TestPoints = cat(5,tau1TestPoints, tau1Min.*A, tau1Max.*A);
  tau2TestPoints = cat(5,tau2TestPoints, tau2, tau2);
  %testPoints{end+1} = cat(5,tau1Min.*A,tau2);
  %testPoints{end+1} = cat(5,tau1Max.*A,tau2);
  
  
end  

%  We will have some duplicates, so let's just keep the unique pairs to
%  test
  
if trim
  tau1Test = zeros(g.N(1),g.N(2),g.N(3),g.N(4),6);
  tau2Test = tau1Test;
  for q = 1:g.N(1)
    for r = 1:g.N(2)
      for s = 1:g.N(3)
        for t = 1:g.N(4)
          test = squeeze([tau1TestPoints(q,r,s,t,:), tau2TestPoints(q,r,s,t,:)])';
          test = unique(test,'rows');
          for z = 1:length(test(:,1))
            tau1Test(q,r,s,t,z)=test(z,1);
            tau2Test(q,r,s,t,z)=test(z,2);
          end
        end
      end
    end
    disp([num2str(q) ' of ' num2str(g.N(1)) ' done'])
  end
  schemeData.tau1test = tau1Test;
  schemeData.tau2test = tau2Test;
else
  schemeData.tau1test = tau1TestPoints;
  schemeData.tau2test = tau2TestPoints;
end

end