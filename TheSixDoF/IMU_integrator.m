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

state = [pos;vel;quat];

end
