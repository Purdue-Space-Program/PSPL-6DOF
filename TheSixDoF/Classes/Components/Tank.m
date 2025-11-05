classdef Tank < RocketComponent

    properties
        Length (1,1) double % Length [m]
        TankDia (1,1) double % Tank Diameter [m]
        FuelOx {mustBeMember(FuelOx,{'Fuel', 'Ox'})} = 'Fuel' % Fuel or Oxidizer ('Fuel', 'Ox')
        LiquidMass (1,1) double % Mass of Liquid Propellant [kg]
    end

    methods
        function obj = Tank(name)
            arguments
                name string
            end
            obj.Name = name;
        end
    end
end