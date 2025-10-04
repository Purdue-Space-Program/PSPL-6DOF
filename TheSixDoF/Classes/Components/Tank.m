classdef Tank < RocketComponent

    properties
    end

    methods
        function obj = Tank(name)
            arguments
                name string
            end

            obj.name = name;
        end
    end

end