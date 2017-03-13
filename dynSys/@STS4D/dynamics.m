function dx = dynamics(obj, ~, x, u, ~)

dx = cell(obj.nx,1);
dims = obj.dims;

returnVector = false;
if ~iscell(x)
  returnVector = true;
  x = num2cell(x);
  u = num2cell(u);
end

for i = 1:length(dims)
  dx{i} = dynamics_cell_helper(obj, x, u, obj.dims, obj.dims(i));
end

if returnVector
  dx = cell2mat(dx);
end
end

function dx = dynamics_cell_helper(obj, x, u, dims, dim)

L1 = obj.L1;
R1 = obj.R1;
R2 = obj.R2;
M1 = obj.M1;
M2 = obj.M2;
g = obj.grav;

%dx{dims==1} = x{dims==2};

%dx{dims==2} = tau1.*(tau1num1./denom1) + tau2.*(tau2num1./denom1) + (num1./denom1);

denom1 = (M1.*R1.^4 + L1.^2.*M2.*R1.^2 - L1.^2.*M2.*R2.^2.*cos(x{1} - x{3}).^2);

num1 = -(M2.*cos(x{1} - x{3}).*sin(x{1} - x{3}).*L1.^2.*R2.^2.*x{2}.^2 + ...
    M2.*sin(x{1} - x{3}).*L1.*R1.^2.*R2.*x{4}.^2 + M2.*g.*sin(x{1}).*L1.*R1.^2 - ...
    M2.*g.*sin(x{3}).*cos(x{1} - x{3}).*L1.*R2.^2 + ...
    M1.*g.*sin(x{1}).*R1.^3);
    
tau1num1 = R1.^2;
tau2num1 =  -cos(x{1} - x{3}).*L1.*R2;
%dx{dims==3} = x{dims==4};

%dx{dims==4} = tau1.*(tau1num2./denom2) + tau2.*(tau2num2./denom2)  + (num2./denom2);
denom2 = M2.*(M1.*R1.^4 + L1.^2.*M2.*R1.^2 - L1.^2.*M2.*R2.^2.*cos(x{1} - x{3}).^2);
num2 =(- L1.^2.*M2.^2.*R2.*g.*sin(x{3}) + L1.^3.*M2.^2.*R2.*x{2}.^2.*sin(x{1} - x{3}) + ...
    L1.^2.*M2.^2.*R2.*g.*sin(x{1}).*cos(x{1} - x{3}) - M1.*M2.*R1.^2.*R2.*g.*sin(x{3}) + ...
    L1.^2.*M2.^2.*R2.^2.*x{4}.^2.*cos(x{1} - x{3}).*sin(x{1} - x{3}) +...
    L1.*M1.*M2.*R1.^2.*R2.*x{2}.^2.*sin(x{1} - x{3}) + L1.*M1.*M2.*R1.*R2.*g.*sin(x{1}).*cos(x{1} - x{3}));
    
tau1num2 = - L1.*M2.*R2.*cos(x{1} - x{3});
tau2num2 = M1.*R1.^2 + L1.^2.*M2;


switch dim
  case 1
    dx = x{dims==2};
  case 2
    dx = u{1}.*(tau1num1./denom1) + u{2}.*(tau2num1./denom1) + (num1./denom1);
    %a = isnan(dx);
    %dx(a) = 0;
  case 3
    dx = x{dims==4};
  case 4
    dx = u{1}.*(tau1num2./denom2) + u{2}.*(tau2num2./denom2) + (num2./denom2);
    %a = isnan(dx);
    %dx(a) = 0;
  otherwise
    error('Only dimension 1-4 are defined for dynamics of STS4D!')
end
end