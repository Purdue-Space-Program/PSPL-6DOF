classdef Sensor
    % The Sensor class defines the generic properties of an imperfect
    % sensor. The sensor class assumes zero-mean gaussian white noise for the
    % sensors to satisfy the assumptions of Kalman filtering. The sensor
    % class assumes digital sensors.
    % The basic Sensor class includes:
    % samplingRate: the rate of data acquisition
    % variance: the variance of the noise in the sensor
    % resolution: the finest percievable resolution of the sensor
    % bias: constant offset in the sensor

    properties
        Name (1,1) string
        SamplingRate (1,1) double = 0.25;
        Variance double
        Resolution (1,1) double 
        Bias (1,1) double = 0
    end

    methods
        function out = Sensor(name, samplingRate, variance, resolution, bias)
            out.Name = name;
            out.SamplingRate = samplingRate;
            out.Variance = variance;
            out.Resolution = resolution;
            out.Bias = bias;
        end

    end
end