classdef Fins < RocketComponent

    properties
        Count (1,1) int8 {mustBeMember(Count, [3, 4])} = 3 % Fin Count

        % add the properties for the airfoil later with mustBeMeber
        
        Airfoil (1,1) string % Airfoil
        Span (1,1) double % Span [m]
        RootChord (1,1) double % Root Chord Length [m]
        TipChord (1,1) double % Tip Chord Length [m]
        Sweep (1,1) double % Sweep [m]
        Thickness (1,1) double % Thickness [m]
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