classdef SolidMotor < RocketComponent
    properties
        Thrust (1,1) double % Thrust [N]
        Length (1,1) double % Length [m]
        Dia (1,1) double % Diameter [m]
        BurnTime (1,1) double % Burn Time [s]
        MassFlow (1,1) double % Mass flow rate [kg/s]
        ExitArea (1,1) double % Exit Area [m²]
        ExitPressure (1,1) double % Exit Pressure [Pa]

    end

    methods
        function motorObj = SolidMotor(name)
            motorObj.Name = name;
        end
    end

end