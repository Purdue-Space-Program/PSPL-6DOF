classdef Rocket
    %Rocket: The rocket class is the class the governs all the
    % properties of the rocket itself. If no inputs are given, the
    % values default to values for CMS.
    %
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

        function rocket = rocket()

            name_prompt = "Input the name of a saved rocket (listed below), '1' to save a new rocket, or '2' to create a temporary object (not saved). ";
            name_struct = dir("TheSixDoF\Inputs\Saved Saved Rockets\*.mat");
            rocket_num = length(name_struct);
            
            input(name_prompt+names,"s")


            rocket.name = name;
            rocket.refArea = refArea;
            rocket.thrust = thrust;
            rocket.exitArea = exitArea;
            rocket.exitPressure = exitPressure;
            rocket.radius = radius;
            rocket.length = length;
            rocket.drag = 
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