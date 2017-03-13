function uOpt = optCtrl(obj, ~, x, deriv, uMode, ~)
% uOpt = optCtrl(obj, t, x, deriv, uMode, dMode, MIEdims)

constraint = 0;


%% Input processing
if nargin < 5
  uMode = 'min';
end

R1 = obj.R1;
R2 = obj.R2;
M1 = obj.M1;
M2 = obj.M2;
L1 = obj.L1;
g = 9.81;


p = deriv;


%% Optimal control
if iscell(deriv)
  %dx{1} = x{2};
  
  %dx{2} = tau1.*(tau1num1./denom1) + tau2.*(tau2num1./denom1) + (num1./denom1);

denom1 = (M1.*R1.^4 + L1.^2.*M2.*R1.^2 - L1.^2.*M2.*R2.^2.*cos(x{1} - x{3}).^2);

num1 = -(M2.*cos(x{1} - x{3}).*sin(x{1} - x{3}).*L1.^2.*R2.^2.*x{2}.^2 + ...
    M2.*sin(x{1} - x{3}).*L1.*R1.^2.*R2.*x{4}.^2 + M2.*g.*sin(x{1}).*L1.*R1.^2 - ...
    M2.*g.*sin(x{3}).*cos(x{1} - x{3}).*L1.*R2.^2 + ...
    M1.*g.*sin(x{1}).*R1.^3);
    
tau1num1 = R1.^2;
tau2num1 =  -cos(x{1} - x{3}).*L1.*R2;
  
  %dx{3} = x{4};
  
  %dx{4} = tau1.*(tau1num2./denom2) + tau2.*(tau2num2./denom2)  + (num2./denom2);
denom2 = M2.*(M1.*R1.^4 + L1.^2.*M2.*R1.^2 - L1.^2.*M2.*R2.^2.*cos(x{1} - x{3}).^2);
num2 =(- L1.^2.*M2.^2.*R2.*g.*sin(x{3}) + L1.^3.*M2.^2.*R2.*x{2}.^2.*sin(x{1} - x{3}) + ...
    L1.^2.*M2.^2.*R2.*g.*sin(x{1}).*cos(x{1} - x{3}) - M1.*M2.*R1.^2.*R2.*g.*sin(x{3}) + ...
    L1.^2.*M2.^2.*R2.^2.*x{4}.^2.*cos(x{1} - x{3}).*sin(x{1} - x{3}) +...
    L1.*M1.*M2.*R1.^2.*R2.*x{2}.^2.*sin(x{1} - x{3}) + L1.*M1.*M2.*R1.*R2.*g.*sin(x{1}).*cos(x{1} - x{3}));
    
tau1num2 = - L1.*M2.*R2.*cos(x{1} - x{3});
tau2num2 = M1.*R1.^2 + L1.^2.*M2;
  
  % tau1.*(p{2}.*(tau1num1./denom1)+p{4}.*(tau1num2./denom2)) +...
  % tau2.*(p{2}.*(tau2num1./denom1)+p{4}.*(tau2num2./denom2)) +...
  % p{1}.*x{2}+p{3}.*x{4} + p{2}.*(num1./denom1) + p{4}.* (num2./denom2)
  
  extraTerms = p{1}.*x{2} + p{3}.*x{4} + p{2}.*num1./denom1+p{4}.*num2./denom2;
  tau1Multiplier = (p{2}.*tau1num1./denom1 + p{4}.*tau1num2./denom2);
  tau2Multiplier = (p{2}.*tau2num1./denom1 + p{4}.*tau2num2./denom2);
  
  

  
  % hamValue = extraTerms + tau1Multiplier.*tau1{1}+tau2Multiplier.*tau2{1};
  % uOpt = cell(obj.nu, 1);
  % uOpt{1} = tau1{1};
  % uOpt{2} = tau2{1};
  
  
  if strcmp(uMode, 'min')
    if constraint == 1
      tau1 = obj.tau1Test;
      tau2 = obj.tau2Test;
      hamValue = 1e6.*ones(size(tau1{1}));
      uOpt = cell(obj.nu, 1);
      uOpt{1} = nan.*ones(size(tau1{1}));
      uOpt{2} = uOpt{1};
      for i = 1:length(tau1)
        hamValueNew = extraTerms + tau1Multiplier.*tau1{i}+tau2Multiplier.*tau2{i};
        update = (hamValueNew<=hamValue);
        uOpt{1}(update) = tau1{i}(update);
        uOpt{2}(update) = tau2{i}(update);
      end
    else
      uOpt{1} = (tau1Multiplier>=0).*obj.T1Min + (tau1Multiplier<0).*obj.T1Max;
      uOpt{2} = (tau2Multiplier>=0).*obj.T2Min + (tau2Multiplier<0).*obj.T2Max;
    end
    
  elseif strcmp(uMode, 'max')
    if constraint == 1
      tau1 = obj.tau1Test;
      tau2 = obj.tau2Test;
      hamValue = -1e6.*ones(size(tau1{1}));
      uOpt = cell(obj.nu, 1);
      uOpt{1} = nan.*ones(size(tau1{1}));
      uOpt{2} = uOpt{1};
      for i = 1:length(tau1)
        hamValueNew = extraTerms + tau1Multiplier.*tau1{i}+tau2Multiplier.*tau2{i};
        update = (hamValueNew>=hamValue);
        uOpt{1}(update) = tau1{i}(update);
        uOpt{2}(update) = tau2{i}(update);
        %       uOpt{1}(update) = update.*tau1{i}+ (1-update).*uOpt{1};
        %       uOpt{2} = update.*tau2{i} + (1-update).*uOpt{2};
      end
    else
      uOpt{1} = (tau1Multiplier>=0).*obj.T1Max + (tau1Multiplier<0).*obj.T1Min;
      uOpt{2} = (tau2Multiplier>=0).*obj.T2Max + (tau2Multiplier<0).*obj.T2Min;
    end
  else
    error('Unknown uMode!')
  end
  
else
    %dx(1) = x(2);
  
  %dx(2) = tau1.*(tau1num1./denom1) + tau2.*(tau2num1./denom1) + (num1./denom1);

denom1 = (M1.*R1.^4 + L1.^2.*M2.*R1.^2 - L1.^2.*M2.*R2.^2.*cos(x(1) - x(3)).^2);

num1 = -(M2.*cos(x(1) - x(3)).*sin(x(1) - x(3)).*L1.^2.*R2.^2.*x(2).^2 + ...
    M2.*sin(x(1) - x(3)).*L1.*R1.^2.*R2.*x(4).^2 + M2.*g.*sin(x(1)).*L1.*R1.^2 - ...
    M2.*g.*sin(x(3)).*cos(x(1) - x(3)).*L1.*R2.^2 + ...
    M1.*g.*sin(x(1)).*R1.^3);
    
tau1num1 = R1.^2;
tau2num1 =  -cos(x(1) - x(3)).*L1.*R2;
  
  %dx(3) = x(4);
  
  %dx(4) = tau1.*(tau1num2./denom2) + tau2.*(tau2num2./denom2)  + (num2./denom2);
denom2 = M2.*(M1.*R1.^4 + L1.^2.*M2.*R1.^2 - L1.^2.*M2.*R2.^2.*cos(x(1) - x(3)).^2);
num2 =(- L1.^2.*M2.^2.*R2.*g.*sin(x(3)) + L1.^3.*M2.^2.*R2.*x(2).^2.*sin(x(1) - x(3)) + ...
    L1.^2.*M2.^2.*R2.*g.*sin(x(1)).*cos(x(1) - x(3)) - M1.*M2.*R1.^2.*R2.*g.*sin(x(3)) + ...
    L1.^2.*M2.^2.*R2.^2.*x(4).^2.*cos(x(1) - x(3)).*sin(x(1) - x(3)) +...
    L1.*M1.*M2.*R1.^2.*R2.*x(2).^2.*sin(x(1) - x(3)) + L1.*M1.*M2.*R1.*R2.*g.*sin(x(1)).*cos(x(1) - x(3)));
    
tau1num2 = - L1.*M2.*R2.*cos(x(1) - x(3));
tau2num2 = M1.*R1.^2 + L1.^2.*M2;
  
  %hamValue = p{1}.*dx{1} + p{2}.*dx{2} + p{3}.*dx{3} + p{4}.*dx{4};
  
  % tau1.*(p{2}.*(tau1num1./denom1)+p{4}.*(tau1num2./denom2)) +...
  % tau2.*(p{2}.*(tau2num1./denom1)+p{4}.*(tau2num2./denom2)) +...
  % p{1}.*x(2)+p{3}.*x(4) + p{2}.*(num1./denom1) + p{4}.* (num2./denom2)
  
  extraTerms = p(1).*x(2) + p(3).*x(4) + p(2).*num1./denom1+p(4).*num2./denom2;
  tau1Multiplier = (p(2).*tau1num1./denom1 + p(4).*tau1num2./denom2);
  tau2Multiplier = (p(2).*tau2num1./denom1 + p(4).*tau2num2./denom2);
  
  
  % hamValue = extraTerms + tau1Multiplier.*tau1{1}+tau2Multiplier.*tau2{1};
  % uOpt = cell(obj.nu, 1);
  % uOpt{1} = tau1{1};
  % uOpt{2} = tau2{1};
  
  if strcmp(uMode, 'min')
    if constraint == 1
      tau1 = obj.tau1Test;
      tau2 = obj.tau2Test;
      uOpt(1) = nan;
      uOpt(2) = nan;
      hamValue = 1e6;
      for i = 1:length(tau1)
        hamValueNew = extraTerms + tau1Multiplier.*tau1(i)+tau2Multiplier.*tau2(i);
        if hamValueNew <= hamValue
          uOpt(1) = tau1(i);
          uOpt(2) = tau2(i);
        end
      end
    else
      uOpt(1) = (tau1Multiplier>=0).*obj.T1Min + (tau1Multiplier<0).*obj.T1Max;
      uOpt(2) = (tau2Multiplier>=0).*obj.T2Min + (tau2Multiplier<0).*obj.T2Max;
    end
    
  elseif strcmp(uMode, 'max')
    if constraint == 1
      tau1 = obj.tau1Test;
      tau2 = obj.tau2Test;
      uOpt(1) = nan;
      uOpt(2) = nan;
      hamValue = -1e6;
      for i = 1:length(tau1)
        hamValueNew = extraTerms + tau1Multiplier.*tau1(i)+tau2Multiplier.*tau2(i);
        if hamValueNew <= hamValue
          uOpt(1) = tau1(i);
          uOpt(2) = tau2(i);
        end
      end
    else
      uOpt(1) = (tau1Multiplier>=0).*obj.T1Max + (tau1Multiplier<0).*obj.T1Min;
      uOpt(2) = (tau2Multiplier>=0).*obj.T2Max + (tau2Multiplier<0).*obj.T2Min;
    end
  else
    error('Unknown uMode!')
  end
end




end