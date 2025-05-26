classdef Rocket
    %Rocket: The rocket class is the class the governs all the
    % properties of the rocket itself. If no inputs are given, the
    % values default to values for CMS.
    %
    % INPUTS:
    % name = the name of the rocket system
    % refArea = reference area for aerodynamics [m^2]
    % thrust = nominal sea level thrust [N]
    % exitArea = exit area of nozzle [m^2]
    % exitPressure = exit pressure [Pa]
    % radius = radius of rocket [m]
    % length = total length of rocket [m]
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
        function rocket = Rocket(name)
            if (nargin == 1)
                rocket.name = name;
            end
        end

        function rocket = CreateRocket(rocket,name,refArea,thrust,exitArea,exitPressure,radius,length)
            arguments
                rocket Rocket
                name (1,1) string
                refArea (1,1) double
                thrust (1,1) double
                exitArea (1,1)
                exitPressure (1,1)
                radius (1,1)
                length (1,1)
            end

            rocket.name = name;
            rocket.refArea = refArea;
            rocket.thrust = thrust;
            rocket.exitArea = exitArea;
            rocket.exitPressure = exitPressure;
            rocket.radius = radius;
            rocket.length = length;
        end

        function drawRocket(rocket)
            [xBody,yBody,zBody] = cylinder(rocket.radius,50);

            figure()
            zBody = zBody * rocket.length;
            surf(xBody,yBody,zBody, 'FaceColor','red', 'LineStyle','none', 'FaceAlpha','1')
            axis equal
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