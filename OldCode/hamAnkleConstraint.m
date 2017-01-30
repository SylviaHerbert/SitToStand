% Computes the optimal Hamiltonian with knee/hip torque input
% while enforcing constraint on ankle torque

grid_dim = size(Ang1);      % size of grid

tau1 = zeros(grid_dim);     % initialize 4D grid of input vals
tau2 = zeros(grid_dim);     % initialize 4D grid of input vals

for i = 1:grid_dim(1)
    for j = 1:grid_dim(2)
        for k = 1:grid_dim(3)
            for l = 1:grid_dim(4)
                ang1 = Ang1(i,j,k,l);   %states
                vel1 = Vel1(i,j,k,l);
                ang2 = Ang2(i,j,k,l);
                vel2 = Vel2(i,j,k,l);

                p1 = p{1}(i,j,k,l);     %costates
                p2 = p{2}(i,j,k,l); 
                p3 = p{3}(i,j,k,l);
                p4 = p{4}(i,j,k,l);

                a1 = R2;
                a2 = - (R2+L1*cos(ang2));
                h1 = grav*(M1*R1*R2 + M2*L1*R2)*sin(ang1) + ...
                    (vel1 + vel2)^2*L1*M2*R2^2*sin(ang2) - ...
                    grav*M2*L1*R2*sin(ang1+ang2)*cos(ang2) + ...
                    M2*L1^2*R2*vel1^2*cos(ang2)*sin(ang2);
                d1 = (L1^2*M2*R2 + M1*R1^2*R2 - ... 
                    L1^2*M2*R2*cos(ang2)^2);

                b1 = -(M2*R2^2 + M2*R2*L1*cos(ang2));
                b2 = -(-M1*R1^2 - M2*R2^2 - M2*L1^2 - 2*M2*R2*L1*cos(ang2));
                h2 = -((M2^2*R2^2*L1*grav + M1*M2*R1*R2^2*grav)*sin(ang1) + ...
                    (-M2^2*R2*L1^2*grav - M1*M2*R1^2*R2*grav)*sin(ang1 + ang2) + ...
                    ((M2^2*R2*L1^3+M1*M2*R1^2*R2*L1)*vel1^2+ ... 
                    (M2^2*R2^3*L1)*(vel1+vel2)^2)*sin(ang2) + ...
                    (M2^2*R2*L1^2*grav + M1*M2*R1*R2*L1*grav)*cos(ang2)*sin(ang1) + ...
                    (M2^2*R2^2*L1^2*(2*vel1^2 + 2*vel1*vel2 + vel2^2))*cos(ang2)*sin(ang2) - ...
                    M2^2*R2^2*L1*grav*sin(ang1 + ang2)*cos(ang2));
                d2 = (M2*R2^2*(M1*R1^2 + M2*L1^2 - M2*L1^2*cos(ang2)^2));

                c1 = (L1.^2.*M2 + M1.*R1.^2 + M2.*R2.^2 + L0.*M2.*R2.*cos(ang1 ...
                    + ang2) + L0.*L1.*M2.*cos(ang1) + L0.*M1.*R1.*cos(ang1) + ...
                    L1.*M2.*R2.*cos(ang2));
                c2 = (M2.*R2.^2 + L0.*M2.*R2.*cos(ang1...
                    + ang2) + L1.*M2.*R2.*cos(ang2));
                e = - M2.*R2.*grav.*sin(ang1 + ang2)...
                    - (L1.*M2.*grav + M1.*R1.*grav).*sin(ang1) + ...
                    - L0.*M2.*R2.*vel1.^2.*sin(ang1 + ang2) - L0.*L1.*M2.*vel1.^2.*sin(ang1) ... 
                    - L0.*M2.*R2.*vel2.^2.*sin(ang1 + ang2) ...
                    - L0.*M1.*R1.*vel1.^2.*sin(ang1) - L1.*M2.*R2.*vel2.^2.*sin(ang2)...  
                    - L0.*M2.*R2.*vel1.*vel2.*sin(ang1 + ang2) ...
                    - L1.*M2.*R2.*vel1.*vel2.*sin(ang2);

                c = [p2*a1/d1 + p4*b1/d2, p2*a2/d1 + p4*b2/d2]';
                d = p1*vel1 + p2*h1/d1 + p3*vel2 + p4*h2/d2;
                b = c1*h1/d1 + c2*h2/d2 + e; 

                A = [a1*c1/d1 + b1*c2/d2, a2*c1/d1 + b2*c2/d2]';

                f = c;
                b1 = tau_a_bound(2) - b;
                b2 = -(tau_a_bound(1) - b);
                
                b1_2 = T0Max - b;
                b2_2 = -(T0Min - b);

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


