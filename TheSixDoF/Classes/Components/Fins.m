classdef Fins < RocketComponent

    properties
        Count (1,1) int8 {mustBeMember(Count, [3, 4])} = 3 % Fin Count
        Airfoil (1,1) string %{mustBeMember(Airfoil, ["Double Wedge", "NACA"])} = "Double Wedge"
        Span (1,1) double % Span [m]
        RootChord (1,1) double % Root Chord Length [m]
        TipChord (1,1) double % Tip Chord Length [m]
        Sweep (1,1) double % Sweep [m]
    end

    methods
        function obj = Fins(name)
            arguments
                name string
            end
            
            obj.Name = name;
        end


    end

end