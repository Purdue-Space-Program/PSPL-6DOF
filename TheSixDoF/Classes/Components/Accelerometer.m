classdef Accelerometer < RocketComponent

    properties
        SamplingRate (1,1) double
        Variance (1, 3) double
        Resolution (1, 1) double
        Bias0 (1, 1) double = 0        % [m/s^2]
        BiasRandomWalkRate = 6.86e-4   % [m/s^2 / sqrt(s)]
        ScaleFactor (1, 1) double = 0
    end

    methods
        function obj = Accelerometer(name)
            arguments
                name string
            end
            
            obj.Name = name;
        end

        function [xyzAccel, tspan] = AccelerometerMeasurement(obj, accel, T,dt)
            arguments
                obj
                accel (:, 3) double
                T double
                dt double = 0;
            end

            bias = obj.Bias0*ones(1,3);
            bias_rw = obj.BiasRandomWalkRate;
            sf = obj.ScaleFactor;

            Var = obj.Variance(1:3);
            Var = reshape(Var, 1, 3);

            n = size(accel, 1);

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

            accelInterp = interp1(tspan, accel, t_imu, 'linear', 'extrap');
            
            accel_imu = accelInterp .* sf + bias + noise;

            xyzAccel = accel_imu;

            % plot:
            figure;
            hold on
            plot(tspan, accel);
            plot(t_imu, accel_imu, 'o', 'MarkerSize', 2);
            legend('$\omega_x$', '$\omega_y$', '$\omega_z$', ...
                '$\omega_x$ gyro', '$\omega_y$ gyro', '$\omega_z$ gyro');
            title("True and Measured Accelerations [m/$s^2$] vs Time (s)")
            xlabel("Time (s)")
            ylabel("Accelerations [rad/$s^2$]")
            hold off
        end
 
    end
end