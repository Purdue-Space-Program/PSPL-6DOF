classdef Tank < RocketComponent

    properties
    end

    methods
        function obj = Tank(name)
            arguments
                name string
            end

            obj.Name = name;
        end
    end

end