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


% Initialize state covariance and measurement covariance
accel_cov = IMU_data.accel_info.Variance;
gyro_cov = IMU_data.gyro_info.Variance';

sig_p = 1; sig_v = 0.01; sig_q = 0.01;
P = diag([sig_p*ones(1,3), sig_v*ones(1,3), sig_q*ones(1,4)].^2);

Q = diag([10,10,10,10,10,10,10,10,10,10]); % Process noise covariance
R = diag(IMU_data.gps_info.Variance);  % Measurement noise covariance


% measurement is only position, and it's world frame
H = eye(3,10);


%% Kalman Filter

% Loop for IMU integration:

output = struct;

Phi = eye(10);

for idx = 1:numel(time)

    % --- State Prediction ---
    % state prediction update:
    curr_time = time(idx);

    [accel,gyro,grav] = ValuesfromIdx(IMU_data,curr_time);

    X(:,:,idx+1) = IMU_Lie_integrator(accel, gyro, grav, X(:,:,idx), dt);

    % Covariance update:

    % ---- Calculate Phi ----
    aw_skew = [ 0,            -a_world(3),  a_world(2);
            a_world(3),    0,           -a_world(1);
           -a_world(2),    a_world(1),   0];

    % Omega skew for the world-frame attitude transition
    w_world = R0 * curr_gyro; % Gyro in world frame
    ww_skew = [ 0,            -w_world(3),  w_world(2);
                w_world(3),    0,           -w_world(1);
               -w_world(2),    w_world(1),   0];
    
    Phi = eye(9);
    Phi(1:3, 1:3) = eye(3) - ww_skew * dt; % Attitude self-transition
    Phi(4:6, 1:3) = -aw_skew * dt;         % Attitude affecting velocity
    Phi(7:9, 4:6) = eye(3) * dt;           % Velocity affecting position

    P = Phi*P*Phi' + Q;

    % --- Measurement Updates ---

    if any(time(idx) == IMU_data.GPS_time)

        gps_idx = find(IMU_data.GPS_time == time(idx));

        z = IMU_data.GPS(gps_idx, :)';

        z_hat = X(1:3, 5);

        % Kalman Gain
        S = H * P * H' + R;
        K = (P * H') / S;

        % update the state
        curr_q = rotm2quat(X(1:3,1:3));
        x_vec = [X(1:3,5); X(1:3,4); curr_q'];

        x_vec = x_vec + K * (z - z_hat);

        % put back in integrator form:
        new_q = x_vec(7:10) / norm(x_vec(7:10));
        X(1:3,1:3,idx+1) = quat2rotm(new_q');
        X(1:3,4,idx+1) = x_vec(4:6); % Velocity
        X(1:3,5,idx+1) = x_vec(1:3); % Position
        
        % Update Covariance
        P = (eye(10) - K * H) * P;
    end 

    % save out the P matrix and the state:
    output.P(:,:,idx) = P;
end

% plot the position over time:

truePosArray = IMU_data.ref_traj;
time = IMU_data.ref_time;
% estimated:
estPosLie = X(1,5,:);
estPosLie = squeeze(estPosLie);

% covariance:
cov = output.P;

cov_xpos = squeeze(cov(1,1,:));


% Rocket Trajectory Plot:
figure;
plot(time,truePosArray(:,1),'r')
hold on
plot(time,estPosLie(1:end-1), 'b')

% sigma bounds:
plot(time,estPosLie(1:end-1)+sqrt(cov_xpos),'r--')
plot(time,estPosLie(1:end-1)-sqrt(cov_xpos),'r--')


xlabel('Dist North (m)');
ylabel('Dist East (m)');
legend('True Trajectory', 'IMU Integration')
