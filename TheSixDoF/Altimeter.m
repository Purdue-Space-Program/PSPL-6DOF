classdef Altimeter < Sensor
    % The Altimeter class is a subclass of Sensor and inherits from the
    % sensor properties. See Sensor for input clarification

    methods
        function alt = Altimeter(name,samplingRate,variance,resolution,bias)
            arguments
                name (1,1) string
                samplingRate (1,1) double
                variance (1,1) double
                resolution (1,1) double = 0
                bias (1,1) double = 0
            end
            alt@Sensor(name,samplingRate,variance,resolution,bias)
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
            % dt = timestep between true height datapoints, must be constant.
            % The timestep must also be smaller than the smallest sensor sampling
            % rate to work properly.
            % vel = velocity input for fluctutations in accuracy based on
            % velocity [m/s]
            arguments
                sensor Sensor
                height double
                dt double = 0 %ignore if no input
                vel double = zeros(1,length(height)) %ignore if no input
            end
            var = sensor.Variance;
            bias = sensor.Bias;

            if (dt == 0)
                alt = zeros(1,length(height));
                for k = 1:length(height)
                    alt(k) = height(k) + bias + randn(1)*sqrt(var) + 0.1*randn(1)*sqrt(var)*sqrt(vel(k));
                end
            else
                sampleSkip = sensor.SamplingRate/dt;
                sampleSkipArray = round(1:sampleSkip:numel(height));
                % initialize to Not a Number (NaN) to ignore other entries
                alt = NaN(1,length(height));

                for k = sampleSkipArray
                    alt(k) = height(k) + bias + randn(1)*sqrt(var) + 0.1*randn(1)*sqrt(var)*sqrt(vel(k));
                end
            end


        end
    end
end