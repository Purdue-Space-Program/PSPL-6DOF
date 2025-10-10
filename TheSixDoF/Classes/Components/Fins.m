classdef Fins < RocketComponent

    properties
        count (1,1) int {mustBeMember(count, [3, 4])}
        airfoil (1,1) string {mustBeMember(airfoil, ["Double Wedge", "NACA"])} = "Double Wedge"
        material (1,1) string
        span (1,1) double
        rootChord (1,1) double
        tipChord (1,1) double
        sweep (1,1) double
    end

    methods
    end

end