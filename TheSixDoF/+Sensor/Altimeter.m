classdef Altimeter < Sensor.Sensor
    % The Altimeter class is a subclass of Sensor and inherits from the
    % sensor properties. See Sensor for input clarification

    methods
        function alt = Altimeter(Name,SamplingRate,Variance,Resolution,Bias,ScaleFactor)
            arguments
                Name (1,1) string
                SamplingRate (1,1) double
                Variance (1,1) double
                Resolution (1,1) double = 0
                Bias (1,1) double = 0
                ScaleFactor (1,1) double = 0
            end
            alt@Sensor.Sensor(Name,SamplingRate,Variance,Resolution,Bias,ScaleFactor)
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
                sensor Sensor.Altimeter
                height double
                dt double = 0 %ignore if no input
                vel double = zeros(1,length(height)) %ignore if no input
            end
            var = sensor.Variance;
            bias = sensor.Bias;
            sf = sensor.ScaleFactor;

            if (dt == 0)
                alt = zeros(1,length(height));
                for k = 1:length(height)
                    alt(k) = height(k) + bias + randn(1)*sqrt(var) + height(k)*sf;
                end
            else
                sampleSkip = sensor.SamplingRate/dt;
                sampleSkipArray = round(1:sampleSkip:numel(height));
                % initialize to Not a Number (NaN) to ignore other entries
                alt = NaN(1,length(height));

                for k = sampleSkipArray
                    alt(k) = height(k) + bias + randn(1)*sqrt(var) + height(k)*sf;
                end
            end
        end
    end
end