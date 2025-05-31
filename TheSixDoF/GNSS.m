classdef GNSS < Sensor
    % The GNSS class is a subclass of Sensor and inherits from the
    % sensor properties. See Sensor for input clarification


    methods

        function gnss = GNSS(name,samplingRate,variance,resolution,bias)
            arguments
                name
                samplingRate (1,1) double
                variance double
                resolution (1,1) double
                bias (1,1) double
            end
            gnss@Sensor(name,samplingRate,variance,resolution,bias)
        end

        function [pos, vel] = GNSSMeasurement(sensor,pos, vel, dt)
            % function to get the altitude measurement from the sensor
            % definition. This can either be run at each timestep, or after
            % the numerical integration
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
                sensor Sensor
                pos (:,3) double
                vel (:,3) double
                dt double = 0 %ignore if no input
            end
            xyVar = sensor.Variance(1);
            zVar = sensor.Variance(2);
            vVar = sensor.Variance(3);

            if (dt == 0)
                alt = zeros(1,length(height));
                for k = 1:length(height)
                    alt(k) = height(k) + randn(1)*sqrt(var) + 0.1*randn(1)*sqrt(var)*sqrt(vel(k));
                end
            else
                sampleSkip = sensor.SamplingRate/dt;
                sampleSkipArray = round(1:sampleSkip:numel(height));
                % initialize to Not a Number (NaN) to ignore other entries
                alt = NaN(1,length(height));

                for k = sampleSkipArray
                    alt(k) = height(k) + randn(1)*sqrt(var) + 0.1*randn(1)*sqrt(var)*sqrt(vel(k));
                end
            end
        end
    end
end