for i = 1:length(data(1,1,1,1,:))
%slices = linspace(0, -pi/2, 30);
%slice = [-pi/4];% pi/2];
%for i = 1:length(slices)
hipAngSlice = 'min';
kneeAngSlice = [-pi/2];
[g3D,data3D]=proj(g,data(:,:,:,:,i),[1 0 0 0],kneeAngSlice); 
%[g2D,data2D]=proj(g3D,data3D,[1 0 0],kneeAngSlice); 

figure(3)
clf

h = visSetIm(g3D, data3D);
%h = visSetIm(g2D, data2D); %, color, level, extraArgs);
xlabel('knee velocity')
ylabel('hip velocity')
%zlabel('hip velocity')
axis([g3D.min(1) g3D.max(1) g3D.min(2) g3D.max(2) g3D.min(3) g3D.max(3)])
pause