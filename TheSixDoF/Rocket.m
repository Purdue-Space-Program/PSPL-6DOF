classdef Rocket
    %Rocket: The rocket class is the superclass the governs all the
    % properties of the rocket simulation. If no inputs are given, the
    % values default to values for CMS.
    %
    % INPUTS:
    % name = the name of the rocket system
    % refArea = reference area for aerodynamics [m^2]
    % thrust = nominal sea level thrust [N]
    % exitArea = exit area of nozzle [m^2]
    % exitPressure = exit pressure [Pa]
    % radius = radius of rocket [m]
    % length = length of rocket [m]
    %
    % OUTPUTS:

    properties
        name (1,1) string = 'CMS'
        refArea (1,1) double = 0.02224
        thrust (1,1) double = 4270.29
        exitArea (1,1) = 0.0070573
        exitPressure (1,1) = 75842.3
        radius (1,1) = 0.0841375
        length (1,1) = 5.123688
    end

    methods
        function obj = Rocket(name)
            if (nargin == 1)
                Rocket.name = name;
                obj = Rocket;
            end
        end

        % function disp(obj)
        %     % function which alters the default behavior of display for
        %     % rocket objects
        % 
        %     % Display rocket name building display string
        %     str = "Rocket name: " + obj.name;
        % 
        %     % Display the string for the rocket
        %     disp(str)
        % end
    end
end