classdef Avionics < RocketComponent
    
    properties
    end

    methods
        function obj = Avionics(name)
            arguments
                name string
            end

            obj.name = name;
        end
    end

end