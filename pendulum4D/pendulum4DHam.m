function hamValue = pendulum4DHam(t, data, deriv, schemeData)
% hamiltonian function for double pendulum
%
% Inputs:
%   schemeData - problem parameters
%     .grid:   grid structure
%     .Tmax:   max torgque
%     R1, R2: radius to center of mass
%     L1, L2: length of each segment
%     M1, M2: mass of each segment
%
% Dynamics:
%   \dot theta1      = omega1
%   \dot omega1      = lots of stuff
%   \dot theta2      = omega2
%   \dot omega2      = lots of stuff
%     |T| <= wMax
%

checkStructureFields(schemeData, 'grid',  'R1', 'R2',...
  'M1','M2','L1','T1Max','T2Max');

g = schemeData.grid;
x{1} = g.xs{1};
x{2} = g.xs{2};
x{3} = g.xs{3};
x{4} = g.xs{4};
p = deriv;

R1 = schemeData.R1;
R2 = schemeData.R2;
L1 = schemeData.L1;
M1 = schemeData.M1;
M2 = schemeData.M2;

tau1test = schemeData.tau1test;
tau2test = schemeData.tau2test;
% T1Max = schemeData.T1Max;
% T1Min = -T1Max;
% T2Max = schemeData.T2Max;
% T2Min = -T2Max;
grav = 9.81;

%dx{1} = x{2};

%dx{2} = tau1.*(tau1num1/denom1) + tau2.*(tau2num1/denom1) + (num1/denom1);
denom1 = (L1.^2.*M2.*R2 + M1.*R1.^2.*R2 - ...
  L1.^2.*M2.*R2.*cos(x{3}).^2);
num1 = grav.*(M1.*R1.*R2 + M2.*L1.*R2).*sin(x{1}) + ...
  (x{2} + x{4}).^2.*L1.*M2.*R2.^2.*sin(x{3}) - ...
  grav.*M2.*L1.*R2.*sin(x{1}+x{3}).*cos(x{3}) + ...
  M2.*L1.^2.*R2.*x{2}.^2.*cos(x{3}).*sin(x{3});
tau1num1 = R2;
tau2num1 = -(R2+L1.*cos(x{3}));

%dx{3} = x{4};

%dx{4} = tau1.*(tau1num2/denom2) + tau2.*(tau2num2/denom2)  + (num2/denom2);
num2 = -((M2.^2.*R2.^2.*L1.*grav + M1.*M2.*R1.*R2.^2.*grav).*sin(x{1}) + ...
  (-M2.^2.*R2.*L1.^2.*grav - M1.*M2.*R1.^2.*R2.*grav).*sin(x{1} + x{3}) + ...
  ((M2.^2.*R2.*L1.^3+M1.*M2.*R1.^2.*R2.*L1).*x{2}.^2+ ...
  (M2.^2.*R2.^3.*L1).*(x{2}+x{4}).^2).*sin(x{3}) + ...
  (M2.^2.*R2.*L1.^2.*grav + M1.*M2.*R1.*R2.*L1.*grav).*cos(x{3}).*sin(x{1}) + ...
  (M2.^2.*R2.^2.*L1.^2.*(2.*x{2}.^2 + 2.*x{2}.*x{4} + x{4}.^2)).*cos(x{3}).*sin(x{3}) - ...
  M2.^2.*R2.^2.*L1.*grav.*sin(x{1} + x{3}).*cos(x{3}));
denom2 = (M2.*R2.^2.*(M1.*R1.^2 + M2.*L1.^2 - M2.*L1.^2.*cos(x{3}).^2));
tau1num2 = -(M2.*R2.^2 + M2.*R2.*L1.*cos(x{3}));
tau2num2= -(-M1.*R1.^2 - M2.*R2.^2 - M2.*L1.^2 - 2.*M2.*R2.*L1.*cos(x{3}));

% p{1}.*dx{1} + p{2}.*dx{2} + p{3}.*dx{3} + p{4}.*dx{4}

% p{1}.*x{2}+p{3}.*x{4} + ...
% p{2}.*(tau1.*(tau1num1/denom1) + tau2.*(tau2num1/denom1)
% +(num1/denom1))+...
% p{4}.*(tau1.*(tau1num2/denom2) + tau2.*(tau2num2/denom2)  +
% (num2/denom2))

% tau1.*(p{2}.*(tau1num1/denom1)+p{4}.*(tau1num2/denom2)) +...
% tau2.*(p{2}.*(tau2num1/denom1)+p{4}.*(tau2num2/denom2)) +...
% p{1}.*x{2}+p{3}.*x{4} + p{2}.*(num1/denom1) + p{4}.* (num2/denom2)



extraTerms = p{1}.*x{2} + p{3}.*x{4} + p{2}.*num1./denom1+p{4}.*num2./denom2;

tau1Multiplier = (p{2}.*tau1num1./denom1 + p{4}.*tau1num2./denom2);

tau2Multiplier = (p{2}.*tau2num1./denom1 + p{4}.*tau2num2./denom2);

  hamValue = extraTerms + tau1Multiplier.*tau1test(:,:,:,:,1)+tau2Multiplier.*tau2test(:,:,:,:,1);
for i = 1:length(tau1test(1,1,1,1,:))
  hamValueNew = extraTerms + tau1Multiplier.*tau1test(:,:,:,:,i)+tau2Multiplier.*tau2test(:,:,:,:,i);
  hamValue = min(hamValue,hamValueNew);
end

hamValue = -hamValue;
%hamValue(isnan(hamValue))=1000;
end