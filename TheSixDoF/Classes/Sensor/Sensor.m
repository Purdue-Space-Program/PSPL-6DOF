classdef Sensor
    % Sensor - The Sensor class defines the generic properties of an
    % imperfect sensor. The sensor class assumes zero-mean gaussian white
    % noise for the sensors to satisfy the assumptions of Kalman filtering.
    %
    % SamplingRate: the time between each data acquisition [s]
    %
    % Variance: the variance of the noise in the sensor. The square of the
    % std. dev. of the sensor noise.
    %
    % Resolution: the finest percievable resolution of the sensor. This
    % corresponds to quantization of the sensor.
    %
    % Bias: constant offset in the sensor outputs from true value.
    %
    % Scale Factor: constant scaling in the sensor outputs from true value.
    % 
    % <a href="https://purdue-space-program.atlassian.net/wiki/spaces/PL/pages/1637023774/Sensor+Class">Documentation</a>

    properties
        Name (1,1) string
        SamplingRate (1,1) double
        Variance double
        Resolution (1,1) double 
        Bias (1,1) double = 0
        ScaleFactor (1,1) double = 0
        
    end

    methods
        function out = Sensor(name, samplingRate, variance, resolution, bias, scaleFactor)
            arguments
                name (1,1) string
                samplingRate (1,1) double   
                variance
                resolution = 0;
                bias = 0;
                scaleFactor = 0;
            end
            out.Name = name;
            out.SamplingRate = samplingRate;
            out.Variance = variance;
            out.Resolution = resolution;
            out.Bias = bias;
            out.ScaleFactor = scaleFactor;
        end
    end
end