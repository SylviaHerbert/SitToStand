classdef STS4D < DynSys
  properties
    R1
    R2
    M1
    M2
    L1
    L0
    grav
    
    T1Max
    T1Min
    T2Max
    T2Min
    TAMax
    TAMin
    
    tau1Test
    tau2Test
    
    dims
  end
  
  methods
    function obj = STS4D(x, R1, R2, M1, M2, L1, L0, grav,...
        T1Max, T1Min, T2Max, T2Min, TAMax, TAMin, dims)
      %
      %
      % Output:
      %   obj       - a DubinsCar object
      
      if numel(x) ~= obj.nx
        error('Initial state does not have right dimension!');
      end
      
      if ~iscolumn(x)
        x = x';
      end
      
      if nargin < 2
        height = 1.72;
        mass = 62;
        M1 = 2*(0.1416*mass);       % mass of thighs
        M2 = (.0694 +.4346)*mass;    % mass of head-arms-trunk
        L0 = .25*height;          % length of segment (shank)
        L1 = .26*height;          % length of segment .4(thigh)
        R1 = .43*L1;          % position of COM along segment (thigh)
        R2 = .6*.4*height;          % position of COM along segment (head-arms-trunk)
      end
      
      if nargin < 9
        T1Max = 107;
        T1Min = -T1Max;
        T2Max = 87;
        T2Min = -60;
        TAMax = 68;
        TAMin = -50;
        grav = 9.81;
      end
      
      if nargin < 15
        dims = 1:4;
      end
      
      % Basic vehicle properties
      obj.pdim = [find(dims == 1) find(dims == 3)]; % Position dimensions
      obj.vdim = [find(dims == 2) find(dims == 4)];   % velocity dimensions
      obj.nx = length(dims);
      obj.nu = 2;
      obj.dims = dims;
      %obj.nd = 3;
      
      obj.x = x;
      obj.xhist = obj.x;
      
%       obj.height = height;
%       obj.mass = mass;
      obj.M1 = M1;       
      obj.M2 = M2;   
      obj.L0 = L0;          
      obj.L1 = L1;          
      obj.R1 = R1;          
      obj.R2 = R2;
      obj.T1Max = T1Max;
      obj.T1Min = T1Min;
      obj.T2Max = T2Max;
      obj.T2Min = T2Min;
      obj.TAMax = TAMax;
      obj.TAMin = TAMin;
      obj.grav = grav;
    end
    
  end % end methods
end % end classdef
