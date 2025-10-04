classdef Output
    % The output class sets all of the outputs of the simulation. This
    % class has options for enabling or disabling any of the plots in
    % simulation output.

    properties
        Figure logical = 1
        Traj logical = 1
        Geoplot logical = 1
    end

    methods (Static)
        function obj = Output(data)
            arguments
                data struct
            end
        end
    end
end