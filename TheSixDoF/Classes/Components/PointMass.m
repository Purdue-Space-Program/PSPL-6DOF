classdef PointMass < RocketComponent

    properties
        Color % Color
    end

    methods
        function obj = PointMass(name)
            arguments
                name string
            end
            obj.Name = name;
        end
    end
end