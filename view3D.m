kneeThetaMin = -2;
kneeThetaMax = 0.23;
figure(1)
hold on
for knee = kneeThetaMin:.05:kneeThetaMax
  clf
        [gOut,dataOut]= proj(g,data,[1 0 0 0],knee);
      
      % get the grid and data all ready for plotting
      [mesh_xs, mesh_data]=gridnd2mesh(gOut,dataOut);
      grid_x = mesh_xs{1};
        grid_y = mesh_xs{2};
        grid_z = mesh_xs{3};
        % create and plot isosurface
        mesh_data = smooth3(mesh_data);
        isurf = isosurface(grid_x,grid_y,grid_z,mesh_data,...
          0);
        h1 = patch(isurf);
          view(-30,30)
          %view(0,0)
  axis vis3d
  axis square
  grid_min = [-1.8091-pi/15, ];
grid_max = [0.0247+pi/15,];
  axis([-0.3289-pi/15,  4.4269+pi/15,...
    0.0013-pi/15, 2.4655+pi/15, -4.9904-pi/15,1.8952+pi/15])
  xlabel('\omega Knee')
  ylabel('\theta Hip')
  zlabel('\omega Hip')
  set(h1, 'FaceColor', 'b', 'EdgeColor', 'none');
  h1.FaceLighting = 'phong';
  title(['\theta Knee = ' num2str(knee)])
  camlight left
  camlight headlight
        drawnow;
        pause
        
end