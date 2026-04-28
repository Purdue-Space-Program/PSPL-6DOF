% this file tests the integration of the IMU data to build a trajectory:
data = load("Rocket_IMU_data.mat");

IMU_data = data.IMU_data; 

time = IMU_data.ref_time;

state = IMU_data.init_cond;

% set up the lie theory:
pos = state(1:3);
vel = state(4:6);
quat = state(7:10);
R = quat2rotm(quat');
X = [R,vel,pos;zeros(2,3),eye(2)];

% integrate the trajectory per time step:
dt = 1/100;

tic
for idx = 1:numel(time)
    curr_time = time(idx);

    [accel,gyro,grav] = ValuesfromIdx(IMU_data,curr_time);

    X(:,:,idx+1) = IMU_Lie_integrator(accel, gyro, grav, X(:,:,idx), dt);

end
toc


% show the integrated trajectory against the truth:

truePosArray = IMU_data.ref_traj;
% estimated:
estPosLie = X(1:3,5,:);
estPosLie = squeeze(estPosLie);
estRotLie = X(1:3,1:3,:);

estQuatLie = rotm2quat(estRotLie);

% Rocket Trajectory Plot:
figure;
plot3(truePosArray(:,1), truePosArray(:,2), truePosArray(:,3), 'g')
hold on
plot3(estPosLie(1,:), estPosLie(2,:), estPosLie(3,:), 'b')
view(43,24);
xlabel('Dist North (m)');
ylabel('Dist East (m)');
zlabel('Height (m)');
legend('True Trajectory', 'IMU Integration')
grid minor;

% compare the quats:
figure;
plot(estQuatLie)
hold on
plot(IMU_data.ref_quat)


figure;
plot(time,truePosArray(:,1),'r')
hold on
plot(time,estPosLie(1,1:end-1), 'b')
plot(IMU_data.GPS_time, IMU_data.GPS(:,1), 'go')


xlabel('Dist North (m)');
ylabel('Dist East (m)');
legend('True Trajectory', 'IMU Integration')
