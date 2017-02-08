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
ang1 = g.xs{1};
vel1 = g.xs{2};
ang2 = g.xs{3};
vel2 = g.xs{4};
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

num1 =(grav.*(M1.*R1.*R2 + M2.*L1.*R2).*sin(ang1) + ...
  (vel1 + vel2).^2.*L1.*M2.*R2.^2.*sin(ang2) - ...
  grav.*M2.*L1.*R2.*sin(ang1+ang2).*cos(ang2) + ...
  M2.*L1.^2.*R2.*vel1.^2.*cos(ang2).*sin(ang2));

denom1 =(L1.^2.*M2.*R2 + M1.*R1.^2.*R2 - ...
  L1.^2.*M2.*R2.*cos(ang2).^2);

num2 = (-(M2.^2.*R2.^2.*L1.*grav + M1.*M2.*R1.*R2.^2.*grav).*sin(ang1)+...
-(-M2.^2.*R2.*L1.^2.*grav - M1.*M2.*R1.^2.*R2.*grav).*sin(ang1 + ang2)+...
-((M2.^2.*R2.*L1.^3+M1.*M2.*R1.^2.*R2.*L1).*vel1.^2+(M2.^2.*R2.^3.*L1).*(vel1+vel2).^2).*sin(ang2)+...
-(M2.^2.*R2.*L1.^2.*grav + M1.*M2.*R1.*R2.*L1.*grav).*cos(ang2).*sin(ang1)+...
-(M2.^2.*R2.^2.*L1.^2.*(2.*vel1.^2 + 2.*vel1.*vel2 + vel2.^2)).*cos(ang2).*sin(ang2)+...
 M2.^2.*R2.^2.*L1.*grav.*sin(ang1 + ang2).*cos(ang2));

denom2 = (M2.*R2.^2.*(M1.*R1.^2 + M2.*L1.^2 - M2.*L1.^2.*cos(ang2).^2));

extraTerms = p{1}.*vel1 + p{3}.*vel2 + p{2}.*num1./denom1+p{4}.*num2./denom2;

tau1Multiplier = (p{2}.*R2./denom1 - p{4}.*(M2.*R2.^2 + M2.*R2.*L1.*cos(ang2))./denom2);

tau2Multiplier = (-p{2}.*(R2+L1.*cos(ang2))./denom1...
    - p{4}.*(-M1.*R1.^2 - M2.*R2.^2 - M2.*L1.^2 - 2.*M2.*R2.*L1.*cos(ang2))./denom2);
  
% hamValue = extraTerms...
%   + (tau1Multiplier>=0).*(tau1Multiplier).*T1Min ...
%   + (tau1Multiplier<0).*(tau1Multiplier).*T1Max ...
%   + (tau2Multiplier>=0).*(tau2Multiplier).*T2Min...
%   + (tau2Multiplier<0).*(tau2Multiplier).*T2Max;

hamValue = extraTerms + tau1Multiplier.*tau1test(:,:,:,:,1)+tau2Multiplier.*tau2test(:,:,:,:,1);
for i = 2:8
  hamValueNew = extraTerms + tau1Multiplier.*tau1test(:,:,:,:,i)+tau2Multiplier.*tau2test(:,:,:,:,i);
  hamValue = min(hamValue,hamValueNew);
end

hamValue = -hamValue;
%hamValue(isnan(hamValue))=1000;
end