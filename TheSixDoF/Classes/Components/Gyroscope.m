classdef Gyroscope < RocketComponent

    properties
        SamplingRate (1,1) double % Airfoil
        Variance double
        Resolution (1, 1) double
        Bias (1, 1) double = 0
        ScaleFactor (1, 1) double = 0
    end

    methods
        function obj = GyroscopeSensor(name)
            arguments
                name string
            end
            
            obj.Name = name;
        end

        function [xyzAngVel] = Gyroscope(obj, angVel, dt)
            arguments
                obj
                angVel (:, 3) double
                dt double = 0;
            end
            bias = obj.Bias;
            sf = obj.ScaleFactor;

            Var = obj.Variance(1:3);

            len = numel(angVel(:, 1));

            omega = angVel(:, 1:3);

            if (dt == 0)
                AngVel = zeros(3, length(omega));
                for k = 1:length(omega)
                    AngVel(k) = omega(k) + bias + randn(3,1) .* sqrt(Var) + omega(k) .* sf;
                end
                xyzAngVel = AngVel;
            else
                sampleSkip = obj.SamplingRate / dt;
                sampleSkipArray = round(1:sampleSkip:len);
            
                AngVel = zeros(3, length(omega));

                for k = sampleSkipArray
                    AngVel(k) = omega(k) + bias + randn(3,1) .* sqrt(Var) + omega(k) .* sf;
                end
                xyzAngVel = AngVel;
            end
        end
    end
end