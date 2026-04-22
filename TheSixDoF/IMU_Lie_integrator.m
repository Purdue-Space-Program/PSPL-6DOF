function state = IMU_Lie_integrator(IMU_data, state, time, dt)
% this function takes the current state and integrates the IMU measurements
% at the given time to produce the updated state. This is the equivalent of
% the prediction step of the Kalman filter equations.
%
% INPUTS:
%
% IMU_data - IMU struct with corrected acceleration, gyro data, etc.
% state - 5x5 matrix of the current state estimate
% time - time at which to start the integration
% dt - integration duration

accel_time = IMU_data.accel_time;
gyro_time = IMU_data.gyro_time;

R0 = state(1:3,1:3);
P0 = state(1:3,4:5);

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
curr_grav  = IMU_data.grav(r_ptr, :)';

% integrate the state
omegaX = gyro(1);
omegaY = gyro(2);
omegaZ = gyro(3);

Omega = [0 -omegaZ omegaY;
        omegaZ 0 -omegaX;
        -omegaY omegaX 0];

theta = norm(gyro);

C1 = (1-theta^2/2-cos(theta)) / (theta^2);
C2 = (theta-sin(theta)) / (theta^3);
C3 = (theta^2/2 - theta^4/24 + cos(theta) - 1) / (theta^4);

AM = zeros(3,2);
AM(3,1) = curr_grav(3);

AN = [accel,zeros(3,1)];

B = [0 1;0 0];

if abs(theta) < 1e-6
    A = Omega*dt;

    sinc_theta = 1 - (theta^2)/6 + (theta^4)/120;
    omc_theta2 = 1/2 - (theta^2)/24 + (theta^4)/720;
    RL = eye(3) + sinc_theta * A + omc_theta2 * A*A;
else
    A = Omega*dt;

    RL = eye(3) + sin(theta)/theta * A + (1-cos(theta))/theta^2* A*A;
end

RR = eye(3);

R = RR*R0*RL;

PM = AM*dt + AM*dt*-B*dt/2;

OmegaT = Omega*dt;
AnT = AN*dt;
BT = B*dt;

PN = AnT + AnT*BT/2 + OmegaT*AnT*(C1*eye(2)+C2*BT)+ ...
    OmegaT*OmegaT*AnT*(C2*eye(2)+C3*BT);


P_mat = RR*R0*PN + (RR*P0+PM)*(eye(2)+BT);

state = [R, P_mat;
    zeros(2,3), eye(2)];
end