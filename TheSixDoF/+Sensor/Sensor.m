classdef Sensor
    % The Sensor class defines the generic properties of an imperfect
    % sensor. The sensor class assumes zero-mean gaussian white noise for the
    % sensors to satisfy the assumptions of Kalman filtering. The sensor
    % class assumes digital sensors.
    % The basic Sensor class includes:
    %
    % samplingRate: the time between each data acquisition [s]
    %
    % variance: the variance of the noise in the sensor. The square of the
    % std. dev. of the sensor noise
    %
    % resolution: the finest percievable resolution of the sensor. This
    % corresponds to quantization of the sensor.
    %
    % bias: constant offset in the sensor. Bias is added to the sensor
    % output

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
            out.Name = name;
            out.SamplingRate = samplingRate;
            out.Variance = variance;
            out.Resolution = resolution;
            out.Bias = bias;
            out.ScaleFactor = scaleFactor;
        end

    end
end