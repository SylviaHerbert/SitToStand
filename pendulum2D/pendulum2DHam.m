function hamValue = pendulum2DHam(t, data, deriv, schemeData)
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

checkStructureFields(schemeData, 'grid','R2','M2','T2Max');

g = schemeData.grid;
ang2 = g.xs{1};
vel2 = g.xs{2};
p = deriv;

R2 = schemeData.R2;
M2 = schemeData.M2;
T2Max = schemeData.T2Max;
T2Min = -T2Max;
grav = 9.81;

hamValue = (p{1}.*vel2)+(p{2}.*(grav./R2).*sin(ang2))...
  +((p{2}.*(1/(M2*R2^2)))>0).*(p{2}.*(1/(M2*R2^2))).*T2Max...
  +((p{2}.*(1/(M2*R2^2)))<=0).*(p{2}.*(1/(M2*R2^2))).*T2Min;


%hamValue = -hamValue;
end