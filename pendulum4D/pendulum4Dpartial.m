function alpha = pendulum4DPartial( ...
  t, data, derivMin, derivMax, schemeData, dim)
% partial function for double pendulum
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

checkStructureFields(schemeData, 'grid',  'R1', 'R2',...
  'M1','M2','L1','T1Max','T2Max')

g = schemeData.grid;
ang1 = g.xs{1};
vel1 = g.xs{2};
ang2 = g.xs{3};
vel2 = g.xs{4};

R1 = schemeData.R1;
R2 = schemeData.R2;
L1 = schemeData.L1;
M1 = schemeData.M1;
M2 = schemeData.M2;

T1Max = schemeData.T1Max;
T2Max = schemeData.T2Max;
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

% extraTerms = p{1}.*vel1 + p{3}.*vel2 + p{2}.*num1./denom1+p{4}.*num2./denom2;
% tau1Multiplier = (p{2}.*R2./denom1 - p{4}.*(M2.*R2.^2 + M2.*R2.*L1.*cos(ang2))./denom2);
% tau2Multiplier = (-p{2}.*(R2+L1.*cos(ang2))./denom1...
%   - p{4}.*(-M1.*R1.^2 - M2.*R2.^2 - M2.*L1.^2 - 2.*M2.*R2.*L1.*cos(ang2))./denom2);

% num1 = (grav.*(M1.*R1.*R2 + M2.*L1.*R2).*sin(ang1) + ...
%     (vel1 + vel2)..^2.*L1.*M2.*R2..^2.*sin(ang2) - ...
%     grav.*M2.*L1.*R2.*sin(ang1+ang2).*cos(ang2) + ...
%     M2.*L1..^2.*R2.*vel1..^2.*cos(ang2).*sin(ang2));
% 
% denom1 = (L1..^2.*M2.*R2 + M1.*R1..^2.*R2 - ... 
%     L1..^2.*M2.*R2.*cos(ang2)..^2);
% 
% num2 = ((M2..^2.*R2..^2.*L1.*grav + M1.*M2.*R1.*R2..^2.*grav).*sin(ang1) + ...
%     (-M2..^2.*R2.*L1..^2.*grav - M1.*M2.*R1..^2.*R2.*grav).*sin(ang1 + ang2) + ...
%     ((M2..^2.*R2.*L1..^3+M1.*M2.*R1..^2.*R2.*L1).*vel1..^2+ ... 
%     (M2..^2.*R2..^3.*L1).*(vel1+vel2)..^2).*sin(ang2) + ...
%     (M2..^2.*R2.*L1..^2.*grav + M1.*M2.*R1.*R2.*L1.*grav).*cos(ang2).*sin(ang1) + ...
%     (M2..^2.*R2..^2.*L1..^2.*(2.*vel1..^2 + 2.*vel1.*vel2 + vel2..^2)).*cos(ang2).*sin(ang2) - ...
%     M2..^2.*R2..^2.*L1.*grav.*sin(ang1 + ang2).*cos(ang2));
%   
% denom2 = (M2.*R2..^2.*(M1.*R1..^2 + M2.*L1..^2 - M2.*L1..^2.*cos(ang2)..^2));

switch dim
%   case 1
%     alpha = abs(vel1);
%     
%   case 2
%     alpha = T1Max.*abs(R2./denom1)...
%       +T2Max.*abs(-(R2+L1.*cos(ang2))./denom1)...
%       +abs(num1./denom1);
%     
%   case 3
%     alpha = abs(vel2);
%   
%   case 4
%     alpha = T1Max.*abs(-(M2.*R2..^2 + M2.*R2.*L1.*cos(ang2))./denom2)...
%       +T2Max.*abs(-(-M1.*R1..^2 - M2.*R2..^2 - M2.*L1..^2 - 2.*M2.*R2.*L1.*cos(ang2))./denom2)...
%       +abs(num2./denom2);
  case 1
    alpha = abs(vel1);
  case 2
    alpha = abs(R2/denom1).*T1Max + (-(R2+L1.*cos(ang2))./denom1).*T2Max + abs(num1./denom1);
  case 3
    alpha = abs(vel2);
  case 4
    alpha = abs(-(M2.*R2.^2 + M2.*R2.*L1.*cos(ang2))./denom2).*T1Max ...
      + abs(-(-M1.*R1.^2 - M2.*R2.^2 - M2.*L1.^2 - 2.*M2.*R2.*L1.*cos(ang2))./denom2).*T2Max...
      + abs(num2./denom2);
end
end