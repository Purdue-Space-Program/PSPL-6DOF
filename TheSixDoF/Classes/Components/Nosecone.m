classdef Nosecone < RocketComponent

    properties
        type (1,1) {mustBeMember(type, ["Von Karman", "Tangent Ogive", "Conical", "Elliptical"])}
        
    end

    methods
    end

end