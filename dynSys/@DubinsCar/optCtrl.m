function uOpt = optCtrl(obj, ~, x, deriv, uMode)
% uOpt = optCtrl(obj, t, y, deriv, uMode)

%% Input processing
if nargin < 5
  uMode = 'min';
end

speed = obj.speed;
p = deriv;
% if ~iscell(deriv)
%   deriv = num2cell(deriv);
% end

%% Optimal control
if iscell(deriv)
  %dx(1) = speed.*cos(x(3));
  %dx(2) = speed.*sin(x(3));
  %dx(3) = u
  uTest = obj.uTest;
uOpt = cell(obj.nu, 1);
uOpt{1} = nan.*ones(size(uTest{1}));

% hamValue = p1.*speed.cos(x(3)) + p2.*speed.sin(x(3)) + p3.*u;
% uOpt = cell(obj.nu, 1);
% uOpt = uTest{1};

  if strcmp(uMode, 'min')
    hamValue = 1e6.*ones(size(uTest{1}));
    for i = 1:length(uTest)
      hamValueNew = p{1}.*speed.*cos(x{3}) + p{2}.*speed.*sin(x{3}) + p{3}.*uTest{i};
      update = (hamValueNew<=hamValue);
      uOpt{1}(update) = uTest{i}(update);
    end
    
  elseif strcmp(uMode, 'max')
    hamValue = -1e6.*ones(size(uTest{1}));
    for i = 1:length(uTest)
      hamValueNew = p{1}.*speed.*cos(x{3}) + p{2}.*speed.*sin(x{3}) + p{3}.*uTest{i};
      update = (hamValueNew>=hamValue);
      uOpt{1}(update) = uTest{i}(update);
    end
  else
    error('Unknown uMode!')
  end
  
else
    %dx(1) = speed.*cos(x(3));
  %dx(2) = speed.*sin(x(3));
  %dx(3) = u
    uTest = obj.uTest;
    hamValue = zeros(1,length(uTest));
    
  for i = 1:length(uTest)
    hamValue(i) = p(1).*speed.*cos(x(3)) + p(2).*speed.*sin(x(3)) + p(3).*uTest(i);
  end
  
  uOpt = zeros(obj.nu, 1);
  if strcmp(uMode, 'max')
    [~,Ind] = max(hamValue(:));
    uOpt(1) = uTest(Ind);
  elseif strcmp(uMode, 'min')
    [~,Ind] = min(hamValue(:));
    uOpt(1) = uTest(Ind);
  else
    error('Unknown uMode!')
  end
end




end