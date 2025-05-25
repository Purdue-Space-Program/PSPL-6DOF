classdef Environment
    % The environment class models the environment system for the rocket,
    % including the sea level height of the location, local gravity, 

    properties
        date datetime = datetime

    end

    methods
        function obj = Environment(inputArg1,inputArg2)
            %UNTITLED5 Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end

        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end