classdef Avionics < RocketComponent
    
    properties
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