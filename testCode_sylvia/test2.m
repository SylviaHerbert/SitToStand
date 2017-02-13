function test2(data,g,tau,nn_intersect,ss_intersect,df_intersect)

%This file makes a gif and an avi file from 3D projections across different
%numbers of psi and v sections. Right now it's set up to do an xyv slice

gif_out_filename = './visualize_STS_BRS_alphaPoint15.gif'; %v_goal.gif';   %Setup filename for gif
avi_out_filename = './visualize_STS_BRS_alphaPoint15.avi'; %v_goal.avi';

% Have to prepare the video writer object first
v = VideoWriter(avi_out_filename);
v.FrameRate = 10;   %Set framerate of playback. 30 is normal.
open(v)   %Opens the file for writing. Make sure to close at the end!!

for i = 1:length(data(1,1,1,1,:))
%slices = linspace(0, -pi/2, 30);
%slice = [-pi/4];% pi/2];
%for i = 1:length(slices)
hipAngSlice = [-pi/2];
kneeAngSlice = [pi/2];
kneeVelSlice = 0;
[g3D,data3D]=proj(g,data(:,:,:,:,i),[1 0 0 0],kneeAngSlice); 
%[g2D,data2D]=proj(g3D,data3D,[1 0 0],kneeAngSlice); 


figure(1)
clf

subplot(1,2,1)
h = visSetIm(g3D, data3D);
%h = visSetIm(g2D, data2D); %, color, level, extraArgs);
xlabel('knee velocity')
ylabel('hip angle')
zlabel('hip velocity')
axis([g3D.min(1) g3D.max(1) g3D.min(2) g3D.max(2) g3D.min(3) g3D.max(3)])
title(['time = ' num2str(tau(i)) ' seconds'])

subplot(1,2,2)
hold on
plot(nn_intersect(:,3), nn_intersect(:,4), 'r*')
plot(df_intersect(:,3), nn_intersect(:,4), 'b*')
plot(ss_intersect(:,3), nn_intersect(:,4), 'g*')
[g2D,data2D]=proj(g3D,data3D,[1 0 0],kneeVelSlice); 
h2 = visSetIm(g2D, data2D);
xlabel('hip angle')
ylabel('hip velocity')
title('slice at transition, alpha = 0.15')

figure(1)
frame = getframe(gcf);   %Get data from figue 1
  image_data = frame2im(frame);   %convert data to image information (this goes straight into avi)
  [imind,cm] = rgb2ind(image_data,256);   %convert image information to index and colormap for gif (don't worry about this part)
  
  if i == 1;   %Make sure this is the loop index == 1
    imwrite(imind,cm,gif_out_filename,'gif', 'Loopcount',inf);   %If you're on the first pass, make a new gif
  else
    imwrite(imind,cm,gif_out_filename,'gif','WriteMode','append');   %If it's not the first pass, append
  end
  
  % Make .avi of same thing
  writeVideo(v,image_data)
end


%Close avi file
 close(v);
end