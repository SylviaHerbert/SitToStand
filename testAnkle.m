function [schemeData] = testAnkle(schemeData)
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

%dx2 = tau1.*(tau1num1/denom1) - tau2.*(tau2num1/denom1) + (num1/denom1);
denom1 = (L1.^2.*M2.*R2 + M1.*R1.^2.*R2 - ...
  L1.^2.*M2.*R2.*cos(x3).^2);
num1 = grav.*(M1.*R1.*R2 + M2.*L1.*R2).*sin(x1) + ...
  (x2 + x4).^2.*L1.*M2.*R2.^2.*sin(x3) - ...
  grav.*M2.*L1.*R2.*sin(x1+x3).*cos(x3) + ...
  M2.*L1.^2.*R2.*x2.^2.*cos(x3).*sin(x3);
tau1num1 = R2;
tau2num1 = (R2+L1.*cos(x3));

%dx3 = x4;

%dx4 = tau1.*(tau1num2/denom2) + tau2.*(tau2num2/denom2)  + (num2/denom2);
num2 = (M2.^2.*R2.^2.*L1.*grav + M1.*M2.*R1.*R2.^2.*grav).*sin(x1) + ...
  (-M2.^2.*R2.*L1.^2.*grav - M1.*M2.*R1.^2.*R2.*grav).*sin(x1 + x3) + ...
  ((M2.^2.*R2.*L1.^3+M1.*M2.*R1.^2.*R2.*L1).*x2.^2+ ...
  (M2.^2.*R2.^3.*L1).*(x2+x4).^2).*sin(x3) + ...
  (M2.^2.*R2.*L1.^2.*grav + M1.*M2.*R1.*R2.*L1.*grav).*cos(x3).*sin(x1) + ...
  (M2.^2.*R2.^2.*L1.^2.*(2.*x2.^2 + 2.*x2.*x4 + x4.^2)).*cos(x3).*sin(x3) - ...
  M2.^2.*R2.^2.*L1.*grav.*sin(x1 + x3).*cos(x3);
denom2 = (M2.*R2.^2.*(M1.*R1.^2 + M2.*L1.^2 - M2.*L1.^2.*cos(x3).^2));
tau1num2 = -(M2.*R2.^2 + M2.*R2.*L1.*cos(x3));
tau2num2= (-M1.*R1.^2 - M2.*R2.^2 - M2.*L1.^2 - 2.*M2.*R2.*L1.*cos(x3));


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


tau1multiplier = (tau1num1./denom1).*tau0num1 + (tau1num2./denom2).*tau0num2;
tau2multiplier = (tau2num1./denom1).*tau0num1 + (tau2num2./denom2).*tau0num2;

% Ankle Constraint that we're solving for:
%  tau0 = tau1.*tau1multiplier - tau2.*tau2multiplier..
%    + (num1/denom1).*tau0num1 + (num2/denom2).*tau0num2 + tau0num3;


%% Solving for test points
%testpoints = which points in control space to test during reachability
%analysis for this set of states

tau1test = cell(8,1);
tau2test = cell(8,1);

for i = 1:2; %for each bound on taus (max&min)
  
  %take that bound for tau1, tau2
  tau1 = tau1Bound(i);
  tau2 = tau2Bound(i);
  
  %plug tau2 bound into ankle constraint, solve for tauAnkleMin + tauAnkleMax
  tau1Test1 = (tauAnkleMin + tau2.*tau2multiplier - (num1./denom1).*tau0num1...
    - (num2./denom2).*tau0num2 - tau0num3)./tau1multiplier;
  tau1Test2 = (tauAnkleMax + tau2.*tau2multiplier - (num1./denom1).*tau0num1...
    - (num2./denom2).*tau0num2 - tau0num3)./tau1multiplier;
  
  %Want to take the intersection between the range of tau1 allowed by the
  %ankle constraint and the range of tau1 itself
  tau1Lim1=max(min(tau1Test1,tau1Test2),tau1Bound(1));
  tau1Lim2=min(max(tau1Test1,tau1Test2),tau1Bound(2));
  
  %%%%% For every grid point want to check if this is an nonempty
  %%%%% intersection between these calculated limits and the given max/min
  %%%%% limits for the input controls
  A = tau1Lim2>=tau1Lim1; %cases when this is a nonempty set
  
  %convert 0's to NaN
  A = double(A);
  A(~A(:))=NaN;
  
  %do the following for indexing purposes
  q=3;
  if i == 1
    q = 0;
  end
  
  %if A is nonempty, make sure we mark the upper and lower bounds for each
  %case
  tau1test{q+i} = tau1Lim1.*A;
  tau2test{q+i}=tau2.*A;
  tau1test{q+i+1}=tau1Lim2.*A;
  tau2test{q+i+1}=tau2.*A;
  
  %do the same thing with tau2
  tau2Test1 = (tau1.*tau1multiplier - tauAnkleMin...
    + (num1./denom1).*tau0num1 + (num2./denom2).*tau0num2 + tau0num3)./tau2multiplier;
  tau2Test2 = (tau1.*tau1multiplier - tauAnkleMax...
    + (num1./denom1).*tau0num1 + (num2./denom2).*tau0num2 + tau0num3)./tau2multiplier;
  
  tau2Lim1=max(min(tau2Test1,tau2Test2),tau2Bound(1));
  tau2Lim2=min(max(tau2Test1,tau2Test2),tau2Bound(2));
  
  B = tau2Lim2>=tau2Lim1;%cases when this is a nonempty set
  %if it is nonempty, make sure we mark the upper and lower bounds for each
  %case
  B = double(B);
  B(~B(:))=NaN; %convert 0's to NaN
  
  tau1test{q+i+2} = tau1.*B;
  tau2test{q+i+2}=tau2Lim1.*B;
  tau1test{q+i+3}=tau1.*B;
  tau2test{q+i+3}=tau2Lim2.*B;
end
schemeData.tau1test = tau1test;
schemeData.tau2test = tau2test;
end