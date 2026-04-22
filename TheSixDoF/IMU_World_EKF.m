% World Frame Kalman Filter for the IMU


%% Setup:

% --- Load the data and set up the initial state ---
data = load("Rocket_IMU_data.mat");

IMU_data = data.IMU_data; 
time = IMU_data.ref_time;
dt = 1/100;
state = IMU_data.init_cond;

% set up the lie group state:
pos = state(1:3);
vel = state(4:6);
quat = state(7:10);
R = quat2rotm(quat');
X = [R,vel,pos;zeros(2,3),eye(2)];

% set up the state and measurement covariance

% get the data from the accelerometer

% Initialize state covariance and measurement covariance
Q = 0.01 * eye(6); % Process noise covariance
R = 4 * eye(3);  % Measurement noise covariance
P = eye(6);        % Initial state covariance

