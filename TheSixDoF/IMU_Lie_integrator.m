function state = IMU_Lie_integrator(accel, gyro, grav, state, dt)
% this function takes the current state and integrates the IMU measurements
% at the given time to produce the updated state. This is the equivalent of
% the prediction step of the Kalman filter equations.
%
% INPUTS:
%
% accel - 3x1 accelometer vector (body frame, specific force)
% gyro - 3x1 gyro vector (body frame)
% state - SE2(3) matrix of current state
% dt - integration duration

R0 = state(1:3,1:3);
P0 = state(1:3,4:5);

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
AM(3,1) = grav(3);

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