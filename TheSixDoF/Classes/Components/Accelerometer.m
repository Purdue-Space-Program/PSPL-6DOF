classdef Accelerometer < RocketComponent

    properties
        SamplingRate (1,1) double % Airfoil
        Variance (1, 3) double
        Resolution (1, 1) double
        Bias (1, 1) double = 0
        ScaleFactor (1, 1) double = 0
    end

    methods
        function obj = Accelerometer(name)
            arguments
                name string
            end
            
            obj.Name = name;
        end

        function [xyzAccel] = AccelerometerMeasurement(obj, accel, dt)
            arguments
                obj Accelerometer
                accel (:, 3) double
                dt double = 0;
            end

            bias = obj.Bias;
            sf = obj.ScaleFactor;

            Var = obj.Variance(1:3);
            Var = reshape(Var, 1, 3);

            n = size(accel, 1);
            AccelerationVector = zeros(n, 3);

            if dt == 0
                sampleIdx = 1:n;
            else
                sampleSkip = obj.SamplingRate / dt;
                sampleIdx = round(1:sampleSkip:n);
            end

            for k = sampleIdx
                noise = randn(1, 3) .* sqrt(Var);
                AccelerationVector(k, :) = accel(k, :) + bias + noise .* sf;
            end

            xyzAccel = AccelerationVector;
        end

        function plotMeasurementHistory(obj, timeArray, accelTrue, accelMeas)
            figure;
            hold on
            plot(timeArray, accelTrue);
            plot(timeArray, accelMeas, 'o', 'MarkerSize', 2);
            legend('X true', 'Y true', 'Z true', 'X sensor', 'Y sensor', 'Z sensor');
            xlabel('Time (s)');
            ylabel('Acceleration on Rocket with Accelerometer Measurement');
            grid on
            hold off
        end
    end
end