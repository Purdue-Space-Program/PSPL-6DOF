classdef Altimeter < Sensor
    % The Altimeter class is a subclass of Sensor and inherits from the
    % sensor properties. See Sensor for input clarification

    methods
        function out = Altimeter(name, samplingRate, variance, resolution, bias)
            out.Name = name;
            out.SamplingRate = samplingRate;
            out.Variance = variance;
            out.Resolution = resolution;
            out.Bias = bias;
        end

        function alt = AltitudeMeasurement(sensor,height, dt, vel)
            % function to get the altitude measurement from the sensor
            % definition. This can either be run at each timestep, or after
            % the numerical integration
            %
            % Required Inputs:
            % sensor = Sensor definition (of type Sensor)
            % height = true height (scalar or vector) [m]
            %
            % Optional Inputs:
            % dt = timestep between height datapoints, must be constant
            % vel = velocity input for fluctutations in accuracy based on
            % velocity [m/s]
            arguments
                sensor Sensor
                height double
                dt double = 0 %ignore if no input
                vel double = 0 %ignore if no input
            end
            var = sensor.Variance;

            if (dt == 0)
                alt = zeros(1,length(height));
                for k = 1:length(height)
                    alt(k) = height(k) + randn(1)*sqrt(var) + randn(1)*sqrt(var)*sqrt(vel(k));
                end
            else
                sampleSkip = sensor.SamplingRate/dt;
                sampleSkipArray = round(1:sampleSkip:numel(height));
                % initialize to Not a Number (NaN) to ignore other entries
                alt = NaN(1,length(height));

                for k = sampleSkipArray
                    alt(k) = height(k) + randn(1)*sqrt(var) + randn(1)*sqrt(var)*sqrt(vel(k));
                end
            end
        end
    end
end