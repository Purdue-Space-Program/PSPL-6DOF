classdef Output
    % The output class sets all of the outputs of the simulation. This
    % class has options for enabling or disabling any of the plots in
    % simulation output.

    properties
        Figure logical = 1
        Traj logical = 1
        Geoplot logical = 1
    end

    methods
        function obj = Output(inputArg1,inputArg2)
            %OUTPUT Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
    end
end