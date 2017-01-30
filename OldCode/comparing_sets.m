T1Max = 1;
T2Max = T1Max;
figure(2)
clf

% gBRS_4D = g;
% dataBRS_4D = data;
% clear g data

%4D BRS from standing, project through knee angle = -pi/2
[gBRS_3D_kneeanglefixed, dataBRS_3D_kneeanglefixed]=...
  proj(gBRS_4D, dataBRS_4D, [1 0 0 0], -pi/2);

%Take result from above, project through all knee velocities to get hip
%info
[gBRS_2D_kneeanglefixed, dataBRS_2D_kneeanglefixed]=...
  proj(gBRS_3D_kneeanglefixed, dataBRS_3D_kneeanglefixed, [1 0 0], 'min');

%Project entire 4D BRS from standing through all knee angles, velocities
[gBRS_2D_all, dataBRS_2D_all] = proj(gBRS_4D, dataBRS_4D, [1 1 0 0], 'min');

%run 2D set forward
[g2Dforward, data2Dforward, data0] = Single_Hips_Forwards(T2Max);
data2Dforward = min(data2Dforward,[],3);

hold on

[~, h0] = contour(g2Dforward.xs{1}, g2Dforward.xs{2}, data2Dforward, [0 0], 'r');
[~, h1] = contour(gBRS_2D_all.xs{1}, gBRS_2D_all.xs{2}, dataBRS_2D_all, [0 0], 'b');
[~, h2] = contour(g2Dforward.xs{1}, g2Dforward.xs{2}, data0, [0 0], 'g');
data0_BRS = shapeRectangleByCorners(gBRS_2D_all, [-pi/80; -.1], [pi/15; .1]);
[~, h3] = contour(gBRS_2D_all.xs{1}, gBRS_2D_all.xs{2}, data0_BRS, [0 0], 'g');
[~, h4] = contour(gBRS_2D_kneeanglefixed.xs{1}, ...
  gBRS_2D_kneeanglefixed.xs{2}, dataBRS_2D_kneeanglefixed, [0 0], 'p');

title(['Projection of Hips for T1Max = T2Max = ' num2str(T1Max) 'Nm, t = 2 sec'],'FontSize',7)
xlabel('\theta Hips','FontSize',7)
ylabel('\omega Hips','FontSize',7)
l=legend('FRS Mode 1','BRS Mode 2', ...
  'Target (Sitting) FRS','Target (Standing) BRS','BRS from Sitting');
