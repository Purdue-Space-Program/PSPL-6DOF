classdef (Abstract) RocketComponent
    properties
        name (1,1) char
        mass (1,1) double {mustBePositive}
        position (1,3) double
        length (1,1) double {mustBePositive}
 
    end

    methods

    end

end