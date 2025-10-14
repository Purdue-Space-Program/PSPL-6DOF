classdef Fins < RocketComponent

    properties
        count (1,1) {mustBeMember(count, [3, 4])} = 3
        airfoil (1,1) string {mustBeMember(airfoil, ["Double Wedge", "NACA"])} = "Double Wedge"
        span (1,1) double
        rootChord (1,1) double
        tipChord (1,1) double
        sweep (1,1) double
    end

    methods
        function obj = Fins(root, tip, span, sweep)
            obj.name = "FinSet";
            obj.rootChord = root;
            obj.tipChord = tip;
            obj.span = span;
            obj.sweep = sweep;

        end


    end

end