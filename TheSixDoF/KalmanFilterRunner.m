data = load("Rocket_IMU_data.mat");

mekf = IMU_ES_MEKF(data);

iekf = IMU_IEKF(data);

mekf_bias = IMU_ES_MEKF_Biases(data);

% plot the results:

% get the covariance for one of them (they are all the same);
cov = mekf.cov_xyz;
time = data.IMU_data.ref_time;

cov_zpos = squeeze(cov(3,3,:));
cov_ypos = squeeze(cov(2,2,:));
cov_xpos = squeeze(cov(1,1,:));



colors = ["#26547c", "#789b73", "#ffd166"];

% Rocket Trajectory Plot:
figure;
sgtitle('EKF Comparison')
subplot(3,1,1)
plot(time,mekf.error(:,3), Color=colors(1), DisplayName='MEKF')
hold on
plot(time,iekf.error(:,3),Color=colors(2), DisplayName='IEKF')
plot(time,mekf_bias.error(:,3),Color=colors(3), DisplayName='MEKF with bias')

% sigma bounds:
plot(time,3*sqrt(cov_zpos),'r--', DisplayName='3-$\sigma$ bounds')
plot(time,-3*sqrt(cov_zpos),'r--', DisplayName='3-$\sigma$ bounds')

xlabel('Time (s)');
ylabel('Up (m)');
legend()
title('Altitude')

subplot(3,1,2)
plot(time,mekf.error(:,2), Color=colors(1), DisplayName='MEKF')
hold on
plot(time,iekf.error(:,2),Color=colors(2), DisplayName='IEKF')
plot(time,mekf_bias.error(:,2),Color=colors(3), DisplayName='MEKF with bias')

% sigma bounds:
plot(time,3*sqrt(cov_ypos),'r--', DisplayName='3-$\sigma$ bounds')
plot(time,-3*sqrt(cov_ypos),'r--', DisplayName='3-$\sigma$ bounds')
%plot(IMU_data.GPS_time, IMU_data.GPS(:,3), 'go')

xlabel('Time (s)');
ylabel('North (m)');
title('East')


subplot(3,1,3)
plot(time,mekf.error(:,1), Color=colors(1), DisplayName='MEKF')
hold on
plot(time,iekf.error(:,1),Color=colors(2), DisplayName='IEKF')
plot(time,mekf_bias.error(:,1),Color=colors(3), DisplayName='MEKF with bias')

% sigma bounds:
plot(time,3*sqrt(cov_xpos),'r--', DisplayName='3-$\sigma$ bounds')
plot(time,-3*sqrt(cov_xpos),'r--', DisplayName='3-$\sigma$ bounds')
%plot(IMU_data.GPS_time, IMU_data.GPS(:,3), 'go')

xlabel('Time (s)');
ylabel('East (m)');
title('North')