classdef (Abstract) RocketComponent
    properties
        Name (1,1) string % Name
        Mass (1,1) double % Mass [kg]
        Position (1,3) double % CoM Pos. [m] (x,y,z)
        Material (1,1) string = "None" % Material
    end

    methods

    end

end