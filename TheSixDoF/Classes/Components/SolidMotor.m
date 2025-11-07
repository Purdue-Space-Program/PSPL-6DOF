classdef SolidMotor < RocketComponent.m
    properties
        Count (1,1) double {mustBePositive} = 1

    end

    methods
        function motorObj = SolidMotor(name)
            motorObj.Name = name;
        end
    end

end