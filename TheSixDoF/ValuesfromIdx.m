function [accel,gyro,grav] = ValuesfromIdx(IMU_data, time)
% this function takes the current state and finds the IMU measurements
% which are closest to that time, replicating a real time integrator.
%
% INPUTS:
%
% IMU_data - IMU struct with corrected acceleration, gyro data, etc.
% time - time at which to start the integration


% find the acceleration and gyro data which are the closest lower value to
% the current time (ZOH):

% Persistent pointers to hold indices between function calls
persistent a_ptr g_ptr r_ptr;

% Initialize pointers if empty or if time has been reset (e.g., new simulation)
if isempty(a_ptr) || time < IMU_data.accel_time(a_ptr)
    a_ptr = 1; g_ptr = 1; r_ptr = 1;
end

% Efficiently increment pointers
% This loop only executes if the time has moved past the current index
while a_ptr < length(IMU_data.accel_time) && IMU_data.accel_time(a_ptr + 1) <= time
    a_ptr = a_ptr + 1;
end

while g_ptr < length(IMU_data.gyro_time) && IMU_data.gyro_time(g_ptr + 1) <= time
    g_ptr = g_ptr + 1;
end

while r_ptr < length(IMU_data.ref_time) && IMU_data.ref_time(r_ptr + 1) <= time
    r_ptr = r_ptr + 1;
end

% Now use the pointers to index the data
accel = IMU_data.accelometer(a_ptr, :)';
gyro  = IMU_data.gyroscope(g_ptr, :)';
grav  = IMU_data.grav(r_ptr, :)';

end