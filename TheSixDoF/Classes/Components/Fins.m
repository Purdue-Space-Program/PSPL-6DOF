classdef Fins < RocketComponent

    properties
        finCount (1,1) int {mustBeMember(finCount, [3, 4])}
        leadingEdgeShape (1,1) string {mustBeMember(leadingEdgeShape, [""])}
        material (1,1) string
        height (1,1) double
        rootChord (1,1) double
        tipChord (1,1) double
    end

    methods
    end

end