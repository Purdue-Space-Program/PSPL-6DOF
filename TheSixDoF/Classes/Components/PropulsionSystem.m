classdef PropulsionSystem < RocketComponent

    properties
        Thrust (1,1) double % Thrust [N]
        BurnTime (1,1) double % Burn Time [s]
        ExitArea (1,1) double % Exit Area [m²]
        ExitPressure (1,1) double % Exit Pressure [Pa]
        Color % Color
    end

    methods
        function obj = PropulsionSystem(name)
            obj.Name = name;
        end
    end
end