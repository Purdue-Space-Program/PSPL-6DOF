classdef Fins < RocketComponent

    properties
        Count (1,1) int8 {mustBeMember(Count, [3, 4])} = 3
        Airfoil (1,1) string {mustBeMember(Airfoil, ["Double Wedge", "NACA"])} = "Double Wedge"
        Span (1,1) double
        RootChord (1,1) double
        TipChord (1,1) double
        Sweep (1,1) double
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