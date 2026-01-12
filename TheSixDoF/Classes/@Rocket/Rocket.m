classdef Rocket < handle

    properties
        Name (1,1) string
        TotalLength % Vehicle Length [m]
        OuterDiameter % Vehicle OD [m]
        NoseLength % Nose Cone Length [m]
        NoseGeometry % Nose Cone Type
        TotalMass % total mass [kg]
        ComponentList dictionary % Component Dictionary
        RASAero_data % RASAero data
        RASAero_data_file_path % The file path of the RASAero data
        CoMOverride % Manual CoM Override
        CoPOverride % Manual CoP Override
    end

    methods
        function obj = Rocket(name)
            obj.Name = name;
            obj.ComponentList = dictionary();
            %obj.CoMOverride = 1.5;
        end


        function addComponent(rocketObj, componentObj)
            arguments
                rocketObj Rocket
                componentObj RocketComponent
            end

            numComponents = numEntries(rocketObj.ComponentList);

            if numComponents == 0
                rocketObj.ComponentList(componentObj.Name) = {componentObj};
            elseif isKey(rocketObj.ComponentList, componentObj.Name)
                %this should be removed for GUI implementation, only added here for testing
                warning("Replacing existing component '%s', do you want to proceed? (Y/N)", componentObj.Name)
                response = input("", "s");
                if strcmp(response, "Y")
                    rocketObj.ComponentList(componentObj.Name) = {componentObj};
                else
                    fprintf("\nStopping . . .\n")
                end
            else
                rocketObj.ComponentList(componentObj.Name) = {componentObj};
            end
        end


        function removeComponent(rocketObj, componentName)
            if isKey(rocketObj.ComponentList, componentName)

                rocketObj.ComponentList = remove(rocketObj.ComponentList, componentName);
            else
                %this should be removed for GUI implementation, only added here for testing
                warning("The component you are trying to remove does not exist, try using a different name")
            end
        end


        function m = getMass(rocketObj)
            arguments
                rocketObj Rocket
            end

            componentArray = values(rocketObj.ComponentList);
            mass = cellfun(@(c) c.mass, componentArray);
            m = sum(mass);
        end


        function RASAero_data = setAeroData(rocket_object)
            arguments
                rocket_object Rocket
            end

            raw_data = readmatrix(rocket_object.RASAero_data_file_path);
            %fprintf("rawData size: %d x %d\n", size(rawData, 1), size(rawData, 2));
            if all(size(raw_data) == [7499, 15])
                RASAero_data = [raw_data(1:300,1:5) raw_data(1:300,8) raw_data(1:300,13:15);
                    raw_data(2501:2800,1:5) raw_data(2501:2800,8) raw_data(2501:2800,13:15);
                    raw_data(5001:5300,1:5) raw_data(5001:5300,8) raw_data(5001:5300,13:15)];
            else
                % data has already been truncated
                RASAero_data = raw_data;
            end

            output_file_path = sprintf("%s_converted_aero_data.csv", rocket_object.Name);
            writematrix(RASAero_data, output_file_path);
        end

        function saveRocket(rocketObj)
            filename = rocketObj.Name;
            filepath = "TheSixDoF" + filesep + "Inputs" + filesep + "Saved Rockets" + filesep + filename + ".mat";
            save(filepath, "rocketObj")
        end


        function A = refArea(rocketObj)
            A = 0;

            components = values(rocketObj.ComponentList);
            for idx = 1:length(components)
                if isa(components(idx), 'Fins')
                    fins = components(idx);
                    A = A + (fins.Thickness * fins.Span * fins.Count);
                end
            end

            A = A + pi * (rocketObj.OuterDiameter / 2)^2;
        end

    end
end