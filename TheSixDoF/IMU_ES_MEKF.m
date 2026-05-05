% Error State Multiplicative Extended Kalman Filter for the IMU

function out = IMU_ES_MEKF(data)

clear ValuesfromIdx

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

sig_p = 2; sig_v = 0.15; sig_theta = 0.5;
P = diag([sig_p*ones(1,3), sig_v*ones(1,3), sig_theta*ones(1,3)].^2);

Q_body = 10*diag([accel_cov,gyro_cov])*dt; % Process noise covariance
R = diag(IMU_data.gps_info.Variance);  % Measurement noise covariance


% measurement is position and velocity, and it's world frame
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
    curr_R = X(1:3, 1:3, idx);

    % Covariance update:

    % Calculate Phi for the given rotation:
    Phi = compute_Phi(curr_R, accel, gyro, dt);

    P = Phi*P*Phi' + Q;
   
    % --- Measurement Updates ---

    if any(time(idx) == IMU_data.GPS_time) && curr_time ~= 0

        gps_idx = find(IMU_data.GPS_time == time(idx));

        z = [IMU_data.GPS(gps_idx, :)';IMU_data.GPS_vel(gps_idx,:)']; 

        z_hat = [X(1:3, 5,idx);X(1:3,4,idx);];

        % store innovations
        innovation = z-z_hat;


        % Kalman Gain
        S = H * P * H' + R;
        K = (P * H') / S;

        NIS(gps_idx) = innovation' * (S \ innovation);

        % get the error state

        delta_xhat = K * (z - z_hat);

        % put the error back into the state
        pos = X(1:3,5,end)+delta_xhat(1:3);
        vel = X(1:3,4,end)+delta_xhat(4:6);


        dq = [1;1/2*delta_xhat(7:9)];
        dq = dq / norm(dq);

        curr_q = rotm2quat(X(1:3,1:3,end));
        updated_q = quatmultiply(curr_q,dq');

        % rebuild the state matrix:
        X(1:3, 5, end) = pos;
        X(1:3, 4, end) = vel;
        X(1:3, 1:3, end) = quat2rotm(updated_q);

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

    out.error = truePosArray-estPosLie(:,1:end-1)';
    out.state = estPosLie(:,1:end-1)';
    out.cov_xyz = squeeze(cov(1:3,1:3,:));


end



% Compute the first order discretization of the phi matrix for the system
function Phi = compute_Phi(R, accel, gyro, dt)

    a_cross = [0 -accel(3) accel(2);
               accel(3) 0 -accel(1);
               -accel(2) accel(1) 0];

    dTheta = (gyro*dt);

    theta_cross = [0 -dTheta(3) dTheta(2);
                   dTheta(3) 0 -dTheta(1);
                   -dTheta(2) dTheta(1) 0];

    F = eye(9);
    F(1:3,4:6) = eye(3)*dt;
    F(4:6,7:9) = -R*a_cross*dt;
    F(7:9,7:9) = eye(3)-theta_cross;

    Phi = F;

end