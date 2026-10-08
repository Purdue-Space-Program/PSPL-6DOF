% Invariant Extended Kalman Filter on SE2(3)

function out = IMU_IEKF(data)

clear ValuesfromIdx

% --- Load the data and set up the initial state ---

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

sig_p = 2; sig_v = 0.15; sig_theta = 0.2;
P = diag([sig_p*ones(1,3), sig_v*ones(1,3), sig_theta*ones(1,3)].^2);

Q_body = 10*diag([accel_cov,gyro_cov])*dt; % Process noise covariance
R = diag(IMU_data.gps_info.Variance);  % Measurement noise covariance

% measurement is position and vel, and it's world frame (identity)
H = eye(6,9);

% convert Q from body into error state frame:
Fi = [zeros(3,6);eye(3) zeros(3);zeros(3) eye(3)];

Q = Fi*Q_body*Fi';


%% Kalman Filter

% Loop for IMU integration:

out = struct;

% initialize the error state:
delta_x = zeros(9,1);

for idx = 1:numel(time)

    % --- State Prediction ---
   
    % state prediction update, find the full state:
    curr_time = time(idx);

    [accel,gyro,grav] = ValuesfromIdx(IMU_data,curr_time);
    X(:,:,idx+1) = IMU_Lie_integrator(accel, gyro, grav, X(:,:,idx), dt);

    % extract the rotation matrix from the full state
    curr_R = X(1:3, 1:3, idx+1);

    % Covariance update:

    % Calculate Phi for the given rotation:
    Phi = compute_Phi_IEKF(accel,gyro,dt);

    P = Phi*P*Phi' + Q;
   
    % --- Measurement Updates ---

    if any(time(idx) == IMU_data.GPS_time) && curr_time ~= 0

        gps_idx = find(IMU_data.GPS_time == time(idx));

        z = IMU_data.GPS(gps_idx, :)';
        z2 = IMU_data.GPS_vel(gps_idx,:)';
     
        z_hat = X(1:3, 5,idx);
        z_hat2 = X(1:3,4,idx);

        

        % get the left invariant innovation
        %Vl = eye(6,9)*(inv(X(:,:,idx+1))*(z-z_hat);

        Vl = [curr_R' * (z-z_hat); curr_R'*(z2-z_hat2)];


        % Kalman Gain
        S = H * P * H' + R;
        K = (P * H') / S;

        % get the error state

        % update the state prediction:

        % wedge the kalman update:
        xi = K*Vl;

        % wedge xi
        Xi = se23_exp(xi);

        X(:,:,idx+1) = X(:,:,idx+1)*Xi;

        % Update Covariance
        P = (eye(9) - K * H) * P*(eye(9)-K*H)' + K*R*K';
    end 

    % save out the P matrix and the state:
    out.P(:,:,idx) = P;
end

% plot the position over time:

truePosArray = IMU_data.ref_traj;
time = IMU_data.ref_time;
% estimated:
estPosLie = X(1:3,5,:);
estPosLie = squeeze(estPosLie);

% covariance:
cov = out.P;

cov_zpos = squeeze(cov(3,3,:));
cov_ypos = squeeze(cov(2,2,:));
cov_xpos = squeeze(cov(1,1,:));

% save the outputs:

out.error = truePosArray-estPosLie(:,1:end-1)';
out.state = estPosLie(:,1:end-1)';
out.cov_xyz = squeeze(cov(1:3,1:3,:));


end


% Compute the exact solution Phi matrix for the system using matrix exp:
function Phi = compute_Phi_IEKF(accel, gyro, dt)
    % (Left-Invariant SE2(3))
    % State Order: [pos (1:3); vel (4:6); rot (7:9)]
    
    O_hat = [ 0,    -gyro(3),  gyro(2);
          gyro(3),  0,    -gyro(1);
         -gyro(2),  gyro(1),  0]; % -omega^ 

    A_hat = [ 0,    -accel(3),  accel(2);
          accel(3),  0,    -accel(1);
         -accel(2),  accel(1),  0]; % -a^
    
    A = zeros(9,9);
    
    % Row 1 (Position dot): -omega^ * p + v
    A(1:3, 1:3) = -O_hat; 
    A(1:3, 4:6) = eye(3);
    
    % Row 2 (Velocity dot): -omega^ * v - a^ * theta
    A(4:6, 4:6) = -O_hat;
    A(4:6, 7:9) = -A_hat;
    
    % Row 3 (Rotation dot): -omega^ * theta
    A(7:9, 7:9) = -O_hat;
    
    % Discretize
    %Phi = expm(A * dt);
    Phi = eye(9)+A*dt;
end


function X = se23_exp(xi)
    % xi is 9x1: [p_x; p_y; p_z; v_x; v_y; v_z; omega_x; omega_y; omega_z;]
    
    phi = xi(7:9); % Rotation vector (axis-angle)
    v   = xi(4:6); % Velocity component
    p   = xi(1:3); % Position component
    
    theta = norm(phi);
    Phi = skew(phi);
    Phi2 = Phi * Phi;
    
    if theta < 1e-6
        % Taylor series for small angles to avoid division by zero
        R = eye(3) + Phi + 0.5 * Phi2;
        V = eye(3) + 0.5 * Phi + (1/6) * Phi2;
    else
        % Rodrigues' Rotation Formula
        R = eye(3) + (sin(theta)/theta) * Phi + ((1-cos(theta))/theta^2) * Phi2;
        
        % Left Jacobian of SO(3)
        V = eye(3) + ((1-cos(theta))/theta^2) * Phi + ((theta-sin(theta))/theta^3) * Phi2;
    end
    
    % Transform velocity and position
    Va = V * v;
    Vb = V * p;
    
    % Construct the 5x5 SE2(3) matrix
    % [ R  Va  Vb ]
    % [ 0   1   0 ]
    % [ 0   0   1 ]
    X = [R,          Va, Vb;
         zeros(1,3), 1,  0;
         zeros(1,3), 0,  1];
end

function s = skew(v)
    % Helper for the 3x3 skew-symmetric matrix
    s = [ 0,    -v(3),  v(2);
          v(3),  0,    -v(1);
         -v(2),  v(1),  0];
end