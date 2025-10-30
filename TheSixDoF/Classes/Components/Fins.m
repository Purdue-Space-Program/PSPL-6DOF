classdef Fins < RocketComponent

    properties
        count (1,1) int8 {mustBeMember(count, [3, 4])} = 3
        airfoil (1,1) string {mustBeMember(airfoil, ["Double Wedge", "NACA"])} = "Double Wedge"
        span (1,1) double
        rootChord (1,1) double
        tipChord (1,1) double
        sweep (1,1) double
    end

    methods
        function obj = Fins(name)
            arguments
                name string
            end
            
            obj.name = name;
        end


    end

end