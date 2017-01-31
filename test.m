load('/Users/sylvia/Documents/MATLAB/SitToStand/Data/standing4DBRS_01_30_2017.mat')

B=data.*(data>-100000);
slices = linspace(0, -pi/2, 30);
%slice = [-pi/4];% pi/2];
for i = 1:length(slices)
[g2D,data2D]=proj(g,B,[1 0 0 0],slices(i));

figure(1)
clf

h = visSetIm(g2D, data2D); %, color, level, extraArgs);
xlabel('knee velocity')
ylabel('hip position')
zlabel('hip velocity')
axis([g2D.min(1) g2D.max(1) g2D.min(2) g2D.max(2) g2D.min(3) g2D.max(3)])
pause
end