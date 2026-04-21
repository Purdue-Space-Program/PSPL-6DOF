% this file tests the integration of the IMU data to build a trajectory:

data = load("Rocket_IMU_data.mat");

IMU_data = data.IMU_data; 

time = IMU_data.ref_time;

state = IMU_data.init_cond;

% integrate the trajectory per time step:

dt = 1/100;

for idx = 1:numel(time)

    curr_time = time(idx);

    state(:,idx+1) = IMU_integrator(IMU_data, state(:,idx), curr_time, dt);

end


% show the integrated trajectory against the truth:

truePosArray = IMU_data.ref_traj;
% estimated:
estPosArray = state(1:3,:);

% Rocket Trajectory Plot:
figure;
plot3(truePosArray(:,1), truePosArray(:,2), truePosArray(:,3), 'g')
hold on
plot3(estPosArray(1,:), estPosArray(2,:), estPosArray(3,:), 'r')
% plot3(posArray(1:endTime / dt,3), posArray(1:endTime / dt,2), zeros(endTime / dt), '--')
% plot3(posArray(1:endTime / dt,3), zeros(endTime / dt), posArray(1:endTime / dt,1), '--')
% plot3(zeros(endTime / dt), posArray(1:endTime / dt,2), posArray(1:endTime / dt,1), '--')
view(43,24);
xlabel('Dist North (m)');
ylabel('Dist East (m)');
zlabel('Height (m)');
legend('True Trajectory', 'IMU Integration')
%axis equal;
grid minor;


function state = IMU_integrator(IMU_data, state, time, dt)
% this function takes the current state and integrates the IMU measurements
% at the given time to produce the updated state. This is the equivalent of
% the prediction step of the Kalman filter equations.
%
% INPUTS:
%
% IMU_data - IMU struct with corrected acceleration, gyro data, etc.
% state - 13x1 vector of the current state estimate
% time - time at which to start the integration
% dt - integration duration

accel = IMU_data.accelometer;
accel_time = IMU_data.accel_time;

gyro = IMU_data.gyroscope;
gyro_time = IMU_data.gyro_time;

% get the current state:
pos = state(1:3);
vel = state(4:6);
quat = state(7:10);

% find the acceleration and gyro data which are the closest lower value to
% the current time (ZOH):

accel_idx = find(accel_time <= time, 1, 'last');
gyro_idx = find(gyro_time <= time, 1, 'last');

curr_accel = accel(accel_idx,:)';
curr_gyro = gyro(gyro_idx,:)';

% find the gravity data which is the closest lower value to
% the current time (ZOH):

grav_idx = find(IMU_data.ref_time <= time, 1, 'last');
curr_grav = IMU_data.grav(grav_idx,:)';


% equations of motion for the update. Quantities are in the body frame, so
% update them first:

R = quat2rotm(quat');
a_world = R*(curr_accel) + curr_grav;

vel = vel + a_world * dt;
pos = pos + vel * dt + 0.5 * a_world * dt^2;

theta = norm(curr_gyro) * dt;

if theta > 1e-8
    axis = curr_gyro / norm(curr_gyro);
    dq = [cos(theta/2);
          axis*sin(theta/2)];
else
    dq = [1; 0; 0; 0];
end

quat = quatmultiply(quat', dq')';
quat = quat / norm(quat);

% omegaX = curr_gyro(1);
% omegaY = curr_gyro(2);
% omegaZ = curr_gyro(3);
% 
% B = [0, -omegaX, -omegaY, -omegaZ;
%      omegaX, 0, omegaZ, -omegaY;
%      omegaY, -omegaZ, 0, omegaX;
%      omegaZ, omegaY, -omegaX, 0];
% 
% quatRates = 0.5 * B * quat;
% 
% quat = quat + quatRates*dt;
% 
state = [pos;vel;quat];

end
