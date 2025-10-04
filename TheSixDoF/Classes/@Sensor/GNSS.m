classdef GNSS < Sensor.Sensor
    % The GNSS class is a subclass of Sensor and inherits from the
    % sensor properties. See Sensor for input clarification


    methods

        function gnss = GNSS(name,samplingRate,variance,resolution,bias,scaleFactor)
            arguments
                name (1,1) string
                samplingRate (1,1) double
                variance (1,3) double
                resolution (1,1) double
                bias (1,1) double
                scaleFactor (1,1) double = 0
            end
            gnss@Sensor.Sensor(name,samplingRate,variance,resolution,bias, scaleFactor)
        end

        function [posOut, velOut] = GNSSMeasurement(sensor,pos, vel, dt)
            % GNSSMeasurement is a method to get the altitude measurement
            % from the sensor definition. This can either be run at each
            % timestep, or after the numerical integration.
            %
            % Required Inputs:
            % sensor = Sensor definition (of type Sensor)
            % pos = true pos (row vector) [m]
            % vel = true vel (row vector) [m/s]
            %
            % Optional Inputs:
            % dt = timestep between GNSS datapoints, must be constant.
            % The timestep must also be smaller than the smallest sensor sampling
            % rate to work properly.
            arguments
                sensor Sensor.GNSS
                pos (:,3) double
                vel (:,3) double
                dt double = 0 %ignore if no input
            end
            xyVar = sensor.Variance(1);
            zVar = sensor.Variance(2);
            vVar = sensor.Variance(3);

            len = numel(pos(:,1));

            if (dt == 0)
                alt = zeros(1,length(height));
                for k = 1:length(height)
                    xyPos(k) = pos(k,1:2) + randn(1)*sqrt(xyVar);
                    zPos(k) = pos(k,3) + randn(1)*sqrt(zVar);
                    velEst(k, :) = vel(k, :) + randn(1, 3) * sqrt(vVar);
                end
            else
                sampleSkip = sensor.SamplingRate/dt;
                sampleSkipArray = round(1:sampleSkip:len);
                % initialize to Not a Number (NaN) to ignore other entries
                xyPos = NaN(len,2);
                zPos = NaN(len, 1);
                velEst = NaN(len, 3);

                for k = sampleSkipArray
                    xyPos(k,:) = pos(k,1:2) + randn(1,2)*sqrt(xyVar);
                    zPos(k,:) = pos(k,3) + randn(1)*sqrt(zVar);
                    velEst(k, :) = vel(k, :) + randn(1, 3) * sqrt(vVar);
                end
                posOut = [xyPos, zPos]; % Combine the 2D and altitude positions
                velOut = velEst; % Output the estimated velocities
            end
        end
    end
end