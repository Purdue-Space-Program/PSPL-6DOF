classdef Gyroscope < RocketComponent

    properties
        SamplingRate (1,1) double % Airfoil
        Variance double
        Resolution (1, 1) double
        Bias (1, 1) double = 0
        ScaleFactor (1, 1) double = 0
    end

    methods
        function obj = Gyroscope(name)
            arguments
                name string
            end
            
            obj.Name = name;
        end

        function [xyzAngVel] = GyroscopeMeasurement(obj, angVel, dt)
            arguments
                obj
                angVel (:, 3) double
                dt double = 0;
            end
            bias = obj.Bias;
            sf = obj.ScaleFactor;

            Var = obj.Variance(1:3);
            Var = reshape(Var, 1, 3);

            n = size(angVel, 1);
            omega = zeros(n, 3);

            if dt == 0
                sampleIdx = 1:n;
            else
                sampleSkip = obj.SamplingRate / dt;
                sampleIdx = round(1:sampleSkip:n);
            end

            for k = sampleIdx
                noise = randn(1, 3) .* sqrt(Var);
                omega(k, :) = angVel(k, :) + bias + noise .* sf;
            end

            xyzAngVel = omega;
        end
        
        function plotMeasurementHistory(obj, timeArray, omegaTrue, omegaMeas)
            % Omega:
            figure;
            hold on
            plot(timeArray, omegaTrue);
            plot(timeArray, omegaMeas, 'o', 'MarkerSize', 2);
            legend('wx', 'wy', 'wz', 'wx gyro', 'wy gyro', 'wz gyro');
            title("True and Measured Angular Velocities [rad/s] vs Time (s)")
            xlabel("Time (s)")
            ylabel("Angular Velocities [rad/s]")
            grid on
            hold off
        end
    end
end