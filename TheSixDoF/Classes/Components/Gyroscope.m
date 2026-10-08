classdef Gyroscope < RocketComponent

    properties
        SamplingRate (1,1) double 
        Variance double 
        Resolution (1, 1) double
        Bias0 (1, 1) double = 0        % [rad/s]
        BiasRandomWalkRate = 1.3e-7    % [rad/s / sqrt(s)]
        ScaleFactor (1, 1) double = 0
    end

    methods
        function obj = Gyroscope(name)
            arguments
                name string
            end
            
            obj.Name = name;
        end

        function [xyzAngVel, tspan] = GyroscopeMeasurement(obj, angVel, T,dt)
            arguments
                obj
                angVel (:, 3) double
                T double
                dt double = 0;
            end

            bias = obj.Bias0*ones(1,3);
            bias_rw = obj.BiasRandomWalkRate;
            sf = obj.ScaleFactor;

            Var = obj.Variance(1:3);
            Var = reshape(Var, 1, 3);

            n = size(angVel, 1);

            % if the sampling rate is zero, match the simulation:
            if obj.SamplingRate == 0
                t_imu = 0:dt:T;
            else
                t_imu = 0:obj.SamplingRate:T;
            end

            % compute the bias random walk:
            N = length(t_imu);

            for k = 2:N
                bias(k,:) = bias(k-1,:) + ...
                    bias_rw .* sqrt(obj.SamplingRate) .*randn(1,3);
            end

            tspan = 0:dt:T+dt;

            % get the gyro data:
            noise = randn(length(t_imu), 3) .* sqrt(Var);

            omegaInterp = interp1(tspan, angVel, t_imu, 'linear', 'extrap');
            
            omega_imu = omegaInterp .* sf + bias + noise;

            xyzAngVel = omega_imu;

            % plot:
            figure;
            hold on
            plot(tspan, angVel);
            plot(t_imu, omega_imu, 'o', 'MarkerSize', 2);
            legend('$\omega_x$', '$\omega_y$', '$\omega_z$', ...
                '$\omega_x$ gyro', '$\omega_y$ gyro', '$\omega_z$ gyro');
            title("True and Measured Angular Velocities [rad/s] vs Time (s)")
            xlabel("Time (s)")
            ylabel("Angular Velocities [rad/s]")
            hold off
        end
        
        function plotMeasurementHistory(obj, timeArray, omegaTrue, omegaMeas)

        end
    end
end