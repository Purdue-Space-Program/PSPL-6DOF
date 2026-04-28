% Invariant Extended Kalman Filter (Right-Invariant) on SE_2(3)

%% Setup

data = load("Rocket_IMU_data.mat");
IMU_data = data.IMU_data;
time     = IMU_data.ref_time;
dt       = 1/100;
state    = IMU_data.init_cond;

pos  = state(1:3);
vel  = state(4:6);
quat = state(7:10);
R0   = quat2rotm(quat');

X_store        = zeros(5,5,numel(time)+1);
X_store(:,:,1) = build_SE2_3(R0, vel(:), pos(:));

% --- Covariances ---
accel_cov = IMU_data.accel_info.Variance(:);
gyro_cov  = IMU_data.gyro_info.Variance(:);

sig_p = 0.5; sig_v = 0.01; sig_theta = 0.01;
P = diag([sig_theta*ones(1,3), sig_v*ones(1,3), sig_p*ones(1,3)].^2);

% Process noise — inflate generously to prevent covariance collapse
B_c = zeros(9,6);
B_c(1:3, 1:3) = eye(3);
B_c(4:6, 4:6) = eye(3);
Q_body = diag([gyro_cov; accel_cov]);
Q = B_c * Q_body * B_c' * dt;

% Add minimum process noise floor to prevent collapse
Q = Q + diag([1e-8*ones(1,3), 1e-6*ones(1,3), 1e-4*ones(1,3)]);

% GPS measurement noise
R_gps = diag(IMU_data.gps_info.Variance) * 2;

%% IEKF Main Loop

output.P = zeros(9,9,numel(time));

for idx = 1:numel(time)

    curr_time = time(idx);
    [accel, gyro, grav] = ValuesfromIdx(IMU_data, curr_time);

    % ----------------------------------------------------------------
    % 1. PROPAGATION
    % ----------------------------------------------------------------

    R_hat = X_store(1:3,1:3,idx);

    X_store(:,:,idx+1) = IMU_Lie_integrator(accel, gyro, grav, X_store(:,:,idx), dt);

    % IEKF state-independent system matrix
    F_c          = zeros(9,9);
    F_c(4:6,1:3) = -hat3(R_hat * accel);
    F_c(7:9,4:6) = eye(3);

    Phi = eye(9) + F_c * dt;
    P   = Phi * P * Phi' + Q;

    % Symmetrize to prevent numerical drift
    P = (P + P') / 2;

    % ----------------------------------------------------------------
    % 2. MEASUREMENT UPDATE (GPS position, world frame)
    % ----------------------------------------------------------------

    gps_idx = find(abs(IMU_data.GPS_time - curr_time) < dt/2, 1);

    if ~isempty(gps_idx) && curr_time > 0

        z_world    = IMU_data.GPS(gps_idx, :)';
        p_hat_curr = X_store(1:3,5,idx+1);
        v_hat_curr = X_store(1:3,4,idx+1);
        R_hat_curr = X_store(1:3,1:3,idx+1);

        % Innovation in world frame
        r = z_world - p_hat_curr;   % 3x1

        % H: position block only — world frame GPS, no attitude coupling
        H = [zeros(3,6), eye(3)];   % 3x9

        S = H * P * H' + R_gps;
        K = (P * H') / S;           % 9x3

        % Full error state
        xi = K * r;                 % 9x1

        delta_phi = xi(1:3);
        delta_v   = xi(4:6);
        delta_p   = xi(7:9);

        % Apply corrections:
        % - Position and velocity: additive (world frame)
        % - Attitude: multiplicative via exp map (on SO(3))
        X_store(1:3,5,idx+1) = p_hat_curr + delta_p;
        X_store(1:3,4,idx+1) = v_hat_curr + delta_v;
        X_store(1:3,1:3,idx+1) = exp_SO3(delta_phi) * R_hat_curr;

        % Joseph-form covariance update
        IKH = eye(9) - K * H;
        P   = IKH * P * IKH' + K * R_gps * K';
        P   = (P + P') / 2;
    end

    output.P(:,:,idx) = P;
end

%% ================================================================
%  Plotting
%% ================================================================

truePosArray = IMU_data.ref_traj;
estPosLie    = squeeze(X_store(1:3,5,:));
cov          = output.P;

cov_xpos = squeeze(cov(7,7,:));
cov_ypos = squeeze(cov(8,8,:));
cov_zpos = squeeze(cov(9,9,:));

figure;
subplot(3,1,1)
plot(time, truePosArray(:,3) - estPosLie(3,1:end-1)')
hold on
plot(time,  3*sqrt(cov_zpos), 'r--')
plot(time, -3*sqrt(cov_zpos), 'r--')
xlabel('Time (s)'); ylabel('(m)');
legend('\delta z', '3\sigma'); title('Altitude Error')

subplot(3,1,2)
plot(time, truePosArray(:,2) - estPosLie(2,1:end-1)')
hold on
plot(time,  3*sqrt(cov_ypos), 'r--')
plot(time, -3*sqrt(cov_ypos), 'r--')
xlabel('Time (s)'); ylabel('(m)');
legend('\delta y', '3\sigma'); title('North Error')

subplot(3,1,3)
plot(time, truePosArray(:,1) - estPosLie(1,1:end-1)')
hold on
plot(time,  3*sqrt(cov_xpos), 'r--')
plot(time, -3*sqrt(cov_xpos), 'r--')
xlabel('Time (s)'); ylabel('(m)');
legend('\delta x', '3\sigma'); title('East Error')

figure;
plot3(truePosArray(:,1), truePosArray(:,2), truePosArray(:,3), 'g')
hold on
plot3(estPosLie(1,1:end-1), estPosLie(2,1:end-1), estPosLie(3,1:end-1), 'b')
view(43,24);
xlabel('Dist North (m)'); ylabel('Dist East (m)'); zlabel('Height (m)');
legend('True Trajectory','IEKF Estimate')
grid minor;

%% ================================================================
%  Helper Functions
%% ================================================================

function X = build_SE2_3(R, v, p)
    X = eye(5);
    X(1:3,1:3) = R;
    X(1:3,4)   = v(:);
    X(1:3,5)   = p(:);
end

function w = hat3(v)
    v = v(:);
    w = [  0,   -v(3),  v(2);
         v(3),    0,   -v(1);
        -v(2),  v(1),    0  ];
end

function R = exp_SO3(phi)
    phi   = phi(:);
    angle = norm(phi);
    if angle < 1e-8
        R = eye(3) + hat3(phi);
    else
        ax = phi / angle;
        K  = hat3(ax);
        R  = eye(3) + sin(angle)*K + (1-cos(angle))*(K*K);
    end
end

function dX = exp_SE2_3(phi, dv, dp)
    phi = phi(:); dv = dv(:); dp = dp(:);
    angle = norm(phi);
    if angle < 1e-8
        dR = eye(3) + hat3(phi);
        V  = eye(3) + 0.5*hat3(phi);
    else
        ax = phi / angle;
        K  = hat3(ax);
        dR = eye(3) + sin(angle)*K + (1-cos(angle))*(K*K);
        V  = eye(3) + ((1-cos(angle))/angle)*K + ((angle-sin(angle))/angle)*(K*K);
    end
    dX = eye(5);
    dX(1:3,1:3) = dR;
    dX(1:3,4)   = V * dv;
    dX(1:3,5)   = V * dp;
end