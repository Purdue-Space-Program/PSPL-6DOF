classdef Accelerometer < RocketComponent

    properties
        SamplingRate (1,1) double % Airfoil
        Variance double
        Resolution (1, 1) double
        Bias (1, 1) double = 0
        ScaleFactor (1, 1) double = 0
    end

    methods
        function obj = AccelerometerComponent(name)
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

            len = numel(accel(:, 1));

            accelVec = accel(1:3);
            if (dt == 0)
                AccelerationVector = zeros(3, length(accelVec));
                for k = 1:length(accelVec)
                    AccelerationVector(k) = accelVec(k) + bias + randn(3,1) .* sqrt(Var) + accelVec(k) .* sf;
                end
                xyzAccel = AccelerationVector;
            else
                sampleSkip = obj.SamplingRate / dt;
                sampleSkipArray = round(1:sampleSkip:len);
            
                AccelerationVector = zeros(3, length(accelVec));

                for k = sampleSkipArray
                    AccelerationVector(k) = accelVec(k) + bias + randn(3,1) .* sqrt(Var) + accelVec(k) .* sf;
                end
                xyzAccel = AccelerationVector;
            end
        end
    end
end