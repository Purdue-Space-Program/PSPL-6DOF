classdef GNSS < Sensor
    % The GNSS class is a subclass of Sensor and inherits from the
    % sensor properties. See Sensor for input clarification

    methods
        function gnss = GNSS(name,samplingRate,variance,resolution,bias,scaleFactor)
            arguments
                name (1,1) string
                samplingRate (1,1) double
                variance (1,6) double
                resolution (1,1) double
                bias (1,1) double
                scaleFactor (1,1) double = 0
            end
            gnss@Sensor(name,samplingRate,variance,resolution,bias, scaleFactor)
        end

        function [posOut, velOut, t_gps] = GNSSMeasurement(sensor,pos,vel, T, dt)
            % GNSSMeasurement is a method to get the position measurement
            % from the sensor definition. This should be run after the
            % numerical integration of the true trajectory
            arguments
                sensor GNSS
                pos (:,3) double
                vel (:,3) double
                T double
                dt double = 0 %ignore if no input
            end

            VarPos = sensor.Variance(1:3);
            VarPos = reshape(VarPos, 1, 3);

            VarVel = sensor.Variance(4:6);
            VarVel = reshape(VarVel, 1, 3);

            % if the sampling rate is zero, match the simulation:
            if sensor.SamplingRate == 0
                t_gps = 0:dt:T;
            else
                t_gps = 0:sensor.SamplingRate:T;
            end

            tspan = 0:dt:T+dt;

            % get the GPS data:
            noisePos = randn(length(t_gps), 3) .* sqrt(VarPos);

            posInterp = interp1(tspan, pos, t_gps, 'linear', 'extrap');
            posOut = posInterp + noisePos;

            noiseVel = randn(length(t_gps), 3) .* sqrt(VarVel);

            velInterp = interp1(tspan, vel, t_gps, 'linear', 'extrap');
            velOut = velInterp + noiseVel;

        end
    end
end
