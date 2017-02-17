function dx = dynamics(obj, ~, x, u, d)
% Dynamics of the Dubins Car
%    \dot{x}_1 = v * cos(x_3)
%    \dot{x}_2 = v * sin(x_3)
%    \dot{x}_3 = w
%   Control: u = w;
%
% Mo Chen, 2016-06-08

if nargin < 5
  d = [0; 0; 0];
end
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

switch dim
  case 1
    dx = obj.speed * cos(x{dims==3});
  case 2
    dx = obj.speed * sin(x{dims==3});
  case 3
    dx = u{1};
    a = isnan(dx);
    dx(a) = 0;
  otherwise
    error('Only dimension 1-3 are defined for dynamics of DubinsCar!')
end
end