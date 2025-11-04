classdef PropulsionSystem < RocketComponent

    properties
        thrust (1,1) double
        burnTime (1,1) double
        fuelMassFlow (1,1) double
        oxMassFlow (1,1) double 
        refArea (1,1) double
        exitPressure (1,1) double
    end

    methods
        function obj = PropulsionSystem(name)
            obj.Name = name;
        end
    end
end