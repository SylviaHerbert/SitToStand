function hamValue = pendulum4DHamAnkle(t, data, deriv, schemeData)
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
L0 = schemeData.L0;

T1Max = schemeData.T1Max;
T1Min = -T1Max;
T2Max = schemeData.T2Max;
T2Min = -T2Max;
grav = 9.81;

tau_a_bound = [-10 10];
T0Max = tau_a_bound(2);
T0Min = tau_a_bound(1);

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
  

% Computes the optimal Hamiltonian with knee/hip torque input
% while enforcing constraint on ankle torque

grid_dim = size(ang1);      % size of grid

tau1 = zeros(grid_dim);     % initialize 4D grid of input vals
tau2 = zeros(grid_dim);     % initialize 4D grid of input vals

for i = 1:grid_dim(1)
    for j = 1:grid_dim(2)
        for k = 1:grid_dim(3)
            for l = 1:grid_dim(4)
                ang1_current = ang1(i,j,k,l);   %states
                vel1_current = vel1(i,j,k,l);
                ang2_current = ang2(i,j,k,l);
                vel2_current = vel2(i,j,k,l);

                p1 = p{1}(i,j,k,l);     %costates
                p2 = p{2}(i,j,k,l); 
                p3 = p{3}(i,j,k,l);
                p4 = p{4}(i,j,k,l);

                a1 = R2;
                a2 = - (R2+L1*cos(ang2_current));
                h1 = grav*(M1*R1*R2 + M2*L1*R2)*sin(ang1_current) + ...
                    (vel1_current + vel2_current)^2*L1*M2*R2^2*sin(ang2_current) - ...
                    grav*M2*L1*R2*sin(ang1_current+ang2_current)*cos(ang2_current) + ...
                    M2*L1^2*R2*vel1_current^2*cos(ang2_current)*sin(ang2_current);
                d1 = (L1^2*M2*R2 + M1*R1^2*R2 - ... 
                    L1^2*M2*R2*cos(ang2_current)^2);

                b1 = -(M2*R2^2 + M2*R2*L1*cos(ang2_current));
                b2 = -(-M1*R1^2 - M2*R2^2 - M2*L1^2 - 2*M2*R2*L1*cos(ang2_current));
                h2 = -((M2^2*R2^2*L1*grav + M1*M2*R1*R2^2*grav)*sin(ang1_current) + ...
                    (-M2^2*R2*L1^2*grav - M1*M2*R1^2*R2*grav)*sin(ang1_current + ang2_current) + ...
                    ((M2^2*R2*L1^3+M1*M2*R1^2*R2*L1)*vel1_current^2+ ... 
                    (M2^2*R2^3*L1)*(vel1_current+vel2_current)^2)*sin(ang2_current) + ...
                    (M2^2*R2*L1^2*grav + M1*M2*R1*R2*L1*grav)*cos(ang2_current)*sin(ang1_current) + ...
                    (M2^2*R2^2*L1^2*(2*vel1_current^2 + 2*vel1_current*vel2_current + vel2_current^2))*cos(ang2_current)*sin(ang2_current) - ...
                    M2^2*R2^2*L1*grav*sin(ang1_current + ang2_current)*cos(ang2_current));
                d2 = (M2*R2^2*(M1*R1^2 + M2*L1^2 - M2*L1^2*cos(ang2_current)^2));

                c1 = (L1.^2.*M2 + M1.*R1.^2 + M2.*R2.^2 + L0.*M2.*R2.*cos(ang1_current ...
                    + ang2_current) + L0.*L1.*M2.*cos(ang1_current) + L0.*M1.*R1.*cos(ang1_current) + ...
                    L1.*M2.*R2.*cos(ang2_current));
                c2 = (M2.*R2.^2 + L0.*M2.*R2.*cos(ang1_current...
                    + ang2_current) + L1.*M2.*R2.*cos(ang2_current));
                e = - M2.*R2.*grav.*sin(ang1_current + ang2_current)...
                    - (L1.*M2.*grav + M1.*R1.*grav).*sin(ang1_current) + ...
                    - L0.*M2.*R2.*vel1_current.^2.*sin(ang1_current + ang2_current) - L0.*L1.*M2.*vel1_current.^2.*sin(ang1_current) ... 
                    - L0.*M2.*R2.*vel2_current.^2.*sin(ang1_current + ang2_current) ...
                    - L0.*M1.*R1.*vel1_current.^2.*sin(ang1_current) - L1.*M2.*R2.*vel2_current.^2.*sin(ang2_current)...  
                    - L0.*M2.*R2.*vel1_current.*vel2_current.*sin(ang1_current + ang2_current) ...
                    - L1.*M2.*R2.*vel1_current.*vel2_current.*sin(ang2_current);

                c = [p2*a1/d1 + p4*b1/d2, p2*a2/d1 + p4*b2/d2]';
                d = p1*vel1_current + p2*h1/d1 + p3*vel2_current + p4*h2/d2;
                b = c1*h1/d1 + c2*h2/d2 + e; 

                A = [a1*c1/d1 + b1*c2/d2, a2*c1/d1 + b2*c2/d2]';

                f = c;
                b1 = T0Max -b;
                b2 = -(T0Min - b);
              

                t = linprog_quiet(f, [A'; -A'], [b1; b2] ,[],[], [T1Min T2Min], [T1Max T2Max]); %solve optimization
                
                tau1(i,j,k,l) = t(1);
                tau2(i,j,k,l) = t(2);
            end
        end
    end
end

hamValue = extraTerms...
  + tau1Multiplier.*tau1 ...
  + tau2Multiplier.*tau2;

hamValue = -hamValue;
end


