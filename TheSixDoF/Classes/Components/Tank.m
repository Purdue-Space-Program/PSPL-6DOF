classdef Tank < RocketComponent

    properties
        Length (1,1) double
        TankDia (1,1) double
        FuelOx {mustBeMember(FuelOx,{'Fuel', 'Prop'})} = 'Fuel'
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