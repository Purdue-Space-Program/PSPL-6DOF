classdef Magnetometer < Sensor.Sensor
    % The Altimeter class is a subclass of Sensor and inherits from the
    % sensor properties. See Sensor for input clarification

    methods
        function alt = Magnetometer(Name,SamplingRate,Variance,Resolution,Bias,ScaleFactor)
            arguments
                Name (1,1) string
                SamplingRate (1,1) double
                Variance (1,3) double
                Resolution (1,1) double = 0
                Bias (1,1) double = 0
                ScaleFactor (1,1) double = 0
            end
            alt@Sensor.Sensor(Name,SamplingRate,Variance,Resolution,Bias,ScaleFactor)
        end

        function [xyzMag] = MagnetometerMeasurement(sensor,env,pos, dt)
            % MagnetometerMeasurement is a method to get the altitude measurement
            % from the sensor definition. This can either be run at each
            % timestep, or after the numerical integration.
            %
            % Required Inputs:
            % sensor = Magnetometer Sensor
            % env = environment 
            % pos = position (lat,long,alt) [m]
            %
            % Optional Inputs: 
            % dt = time between Magnetometer datapoints,
            % must be constant. The simulation timestep must also be
            % smaller than the smallest sensor sampling rate to work
            % properly.
            arguments
                sensor Sensor.Magnetometer
                env Env.Environment
                pos (:,3) double
                dt double = 0 %ignore if no input
            end
            xyVar = sensor.Variance(1);
            zVar = sensor.Variance(2);
            vVar = sensor.Variance(3);

            len = numel(pos(:,1));

            height = pos(:,3);
            lat = pos(:,1);
            long = pos(:,2);

            decimalYear = decyear(env.Date);

            if (dt == 0)
                alt = zeros(1,length(height));
                for k = 1:length(height)
                    xyzMag = wrldmagm(height(k),lat(k),long(k),decimalYear);
                end
            else
                sampleSkip = sensor.SamplingRate/dt;
                sampleSkipArray = round(1:sampleSkip:len);
                % initialize to Not a Number (NaN) to ignore other entries
                xyzMag = NaN(len,3);

                for k = sampleSkipArray
                    xyzMag = wrldmagm(height(k),lat(k),long(k),decimalYear);
                end
            end
        end
    end
end