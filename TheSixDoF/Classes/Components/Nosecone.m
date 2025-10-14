classdef Nosecone < RocketComponent

    properties
        type (1,1) {mustBeMember(type, ["Von Karman", "Tangent Ogive", "Conical", "Elliptical"])} = "Von Karman"
        
    end

    methods
        function obj = Nosecone(name)
            obj.name = name;
        end
    end

end