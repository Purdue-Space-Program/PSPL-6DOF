classdef Gyroscope < Sensor
    % The Gyroscope class is a subclass of Sensor and inherits from the
    % sensor properties. See Sensor for input clarification

    methods
        function alt = Gyroscope(Name,SamplingRate,Variance,Resolution,Bias,ScaleFactor)
            arguments
                Name (1,1) string
                SamplingRate (1,1) double
                Variance (1,3) double
                Resolution (1,1) double = 0
                Bias (1,1) double = 0
                ScaleFactor (1,1) double = 0
            end
            alt@Sensor(Name,SamplingRate,Variance,Resolution,Bias,ScaleFactor)
        end

        function [xyzAngVel] = GyroscopeMeasurement(sensor,angVel, dt)
            % GyroscopeMeasurement is a method to get the angular velocity measurement
            % from the sensor definition. This can either be run at each
            % timestep, or after the numerical integration.
            %
            % Required Inputs:
            % sensor = Gyroscope Sensor
            % angVel = angular velocity (x,y,z) [rad/s]
            %
            % Optional Inputs: 
            % dt = time between Gyroscope datapoints,
            % must be constant. The simulation timestep must also be
            % smaller than the smallest sensor sampling rate to work
            % properly.
            arguments
                sensor Sensor.Gyroscope
                angVel (:,3) double
                dt double = 0 %ignore if no input
            end
            bias = sensor.Bias;
            sf = sensor.ScaleFactor;

            xVar = sensor.Variance(1);
            yVar = sensor.Variance(2);
            zVar = sensor.Variance(3);

            len = numel(angVel(:,1));

            omegaX = angVel(:,1);
            omegaY = angVel(:,2);
            omegaZ = angVel(:,3);

            if (dt == 0)
                xAngVel = zeros(1,length(omegaX));
                yAngVel = zeros(1,length(omegaX));
                zAngVel = zeros(1,length(omegaX));
                for k = 1:length(omegaX)
                    xAngVel(k) = omegaX(k) + bias + randn(1)*sqrt(xVar) + omegaX(k)*sf;
                    yAngVel(k) = omegaY(k) + bias + randn(1)*sqrt(yVar) + omegaY(k)*sf;
                    zAngVel(k) = omegaZ(k) + bias + randn(1)*sqrt(zVar) + omegaZ(k)*sf;
                end
                xyzAngVel=[xAngVel;yAngVel;zAngVel];
            else
                sampleSkip = sensor.SamplingRate/dt;
                sampleSkipArray = round(1:sampleSkip:len);
                % initialize to Not a Number (NaN) to ignore other entries
                xAngVel = zeros(1,length(omegaX));
                yAngVel = zeros(1,length(omegaX));
                zAngVel = zeros(1,length(omegaX));

                for k = sampleSkipArray
                    xAngVel(k) = omegaX(k) + bias + randn(1)*sqrt(xVar) + omegaX(k)*sf;
                    yAngVel(k) = omegaY(k) + bias + randn(1)*sqrt(yVar) + omegaY(k)*sf;
                    zAngVel(k) = omegaZ(k) + bias + randn(1)*sqrt(zVar) + omegaZ(k)*sf;
                end
                xyzAngVel=[xAngVel;yAngVel;zAngVel];
            end
        end
    end
end