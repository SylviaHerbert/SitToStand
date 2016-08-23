%get 2D version of Mode 2 at initial time step with knees still and locked in place
h2 = figure(2);
camlight left
for i = -.03:.0001:0
  %for i = -2.3091:.1:.5247
  %for i = -3.789:.1:4.9269
  clf
  [g3D,data3D]=proj(gMode2,dataMode2(:,:,:,:,end),[0 1 0 0],[i]);
  h2 = visSetIm(g3D, data3D,'g',0,0,0);
  disp(i)
  drawnow;
  pause
end


% [gMode2_comp_test, dataMode2_comp_test] = proj(gMode2, dataMode2(:,:,:,:,end), [1 1 0 0], [-pi/2-.1,-.1;-pi/2+.1,.1]);
% 
% %take final time step for Mode 1
% dataMode1_comp = dataMode1(:,:,end);

% set matrix that is nonzeros 
% dataMode1_comp_inSet = double(dataMode1_comp<=0);
% dataMode2_comp_inSet = double(dataMode2_comp<=0);
% 
% Intersect = double(dataMode1_comp<=0 & dataMode2_comp<=0)