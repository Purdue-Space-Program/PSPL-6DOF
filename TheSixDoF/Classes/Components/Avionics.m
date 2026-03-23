classdef Avionics < RocketComponent
    
    properties
        Magnetometer (1, 1) Magnetometer
        Gyroscope (1, 1) Gyroscope
        GNSS (1, 1) GNSS
        Altimeter (1, 1) Altimeter
    end

    methods
        function obj = Avionics(name)
            arguments
                name string
            end

            obj.Name = name;
        end
    end

end