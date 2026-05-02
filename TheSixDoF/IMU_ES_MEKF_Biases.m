    % Error State Multiplicative Extended Kalman Filter with Bias estimation

    function out = IMU_ES_MEKF_Biases(data)
    
    
    %% Setup:

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
    accel_rw = IMU_data.accel_info.BiasRandomWalkRate;
    gyro_rw = IMU_data.gyro_info.BiasRandomWalkRate;
    
    
    sig_p = 2; sig_v = 0.01; sig_theta = 0.2; sig_accel = 0.1; sig_gyro = 0.1;
    P = diag([sig_p*ones(1,3), sig_v*ones(1,3), sig_theta*ones(1,3),...
        sig_accel*ones(1,3),sig_gyro*ones(1,3)].^2);
    
    Q_body = 10*diag([accel_cov,gyro_cov,ones(1,3)*accel_rw,...
        ones(1,3)*gyro_rw])*dt; % Process noise covariance

    R = diag(IMU_data.gps_info.Variance);  % Measurement noise covariance
    
    
    % measurement is position and velocity, and it's world frame
    H = eye(6,15);
    
    % convert Q from body into error state frame:
    Fi = [zeros(3,12);eye(3) zeros(3,9);zeros(3) eye(3) zeros(3,6); ...
        zeros(3,6) eye(3) zeros(3); zeros(3,9) eye(3)];
    
    Q = Fi*Q_body*Fi';
    
    
    
    %% Kalman Filter
    
    % Loop for IMU integration:
    
    out = struct;
    
    % initialize the error state:
    
    % error state: pos, vel, orientation, accel bias, gyro bias
    delta_x = zeros(15,1);
    
    accel_bias = 0;
    gyro_bias = 0;
    
    for idx = 1:numel(time)
    
        % --- State Prediction ---
       
        % state prediction update, find the full state:
        curr_time = time(idx);
    
        [accel,gyro,grav] = ValuesfromIdx(IMU_data,curr_time);
    
        % subtract the estimated accel bias and gyro bias from measurements:
        accel = accel - accel_bias;
        gyro = gyro - gyro_bias;
    
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
            innovations(gps_idx, :) = z-z_hat;
    
    
            % Kalman Gain
            S = H * P * H' + R;
            K = (P * H') / S;
    
            innov_bounds(gps_idx, :) = 3 * sqrt(diag(S))';
    
            % get the error state
    
            delta_xhat = K * (z - z_hat);
    
            % get the accel bias and gyro bias:
            accel_bias = accel_bias + delta_xhat(10:12);
            gyro_bias = gyro_bias + delta_xhat(13:15);
    
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
            P = (eye(15) - K * H) * P*(eye(15)-K*H)' + K*R*K';
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
    
    cov_accel = squeeze(cov(10,10,:));
    cov_gyro = squeeze(cov(13,13,:));

    out.error = truePosArray-estPosLie(:,1:end-1)';
    out.state = estPosLie(:,1:end-1)';
    out.cov_xyz = squeeze(cov(1:3,1:3,:));
    
    
    % Rocket Trajectory Plot:
    % figure;
    % plot3(truePosArray(:,1), truePosArray(:,2), truePosArray(:,3), 'g')
    % hold on
    % plot3(estPosLie(1,:), estPosLie(2,:), estPosLie(3,:), 'b')
    % view(43,24);
    % xlabel('Dist North (m)');
    % ylabel('Dist East (m)');
    % zlabel('Height (m)');
    % legend('True Trajectory', 'IMU Integration')
    % grid minor;
    
    % Innovations:
    % figure;
    % subplot(2,1,1);
    % plot(IMU_data.GPS_time, innovations(:, 3), 'b'); % Altitude Innovation
    % hold on;
    % plot(IMU_data.GPS_time, innov_bounds(:, 3), 'r--');
    % plot(IMU_data.GPS_time, -innov_bounds(:, 3), 'r--');
    % title('Altitude Innovation ($z - z_{hat}$)');
    % ylabel('Error (m)');
    % grid on;
    % 
    % subplot(2,1,2);
    % plot(IMU_data.GPS_time, innovations(:, 6), 'b'); % Vertical Velocity Innovation
    % hold on;
    % plot(IMU_data.GPS_time, innov_bounds(:, 6), 'r--');
    % plot(IMU_data.GPS_time, -innov_bounds(:, 6), 'r--');
    % title('Vertical Velocity Innovation');
    % ylabel('Error (m/s)');
    % xlabel('Time (s)');
    % grid on;
    
    
    % % compare the quats:
    % estRotLie = X(1:3,1:3,:);
    % estQuatLie = rotm2quat(estRotLie);
    % 
    % figure;
    % plot(estQuatLie, 'DisplayName', 'Estimated Quat')
    % hold on
    % plot(IMU_data.ref_quat, 'DisplayName', 'True Quat')
    % legend();

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
    
        F = eye(15);
        F(1:3,4:6) = eye(3)*dt;
        F(4:6,7:9) = -R*a_cross*dt;
        F(7:9,7:9) = eye(3)-theta_cross;
        F(4:6,10:12) = -R*dt;
        F(7:9,13:15) = -eye(3)*dt;
    
        Phi = F;
    
    end