classdef (Abstract) RocketComponent
    properties
        Name (1,1) string % Name
        Mass (1,1) double % Mass [kg]
        Position (1,1) double % Position [m] (x,y,z)
        Color
        Material (1,1) string = "None"
 
    end

    methods

    end

end