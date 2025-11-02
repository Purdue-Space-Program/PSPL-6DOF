classdef Nosecone < RocketComponent

    properties
        Type (1,1) {mustBeMember(Type, ["Von Karman", "Tangent Ogive", "Conical", "Elliptical"])} = "Von Karman"
        
    end

    methods
        function obj = Nosecone(name)
            arguments
               name string
            end
            
            obj.name = name;
        end
    end

end