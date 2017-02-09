% plot Intersection

%load intersection

figure(1)
hold on
plot(nn_intersect(:,3), nn_intersect(:,4), 'r*')
plot(df_intersect(:,3), nn_intersect(:,4), 'b*')
plot(ss_intersect(:,3), nn_intersect(:,4), 'g*')

xlabel('theta_h')
ylabel('omega_h')