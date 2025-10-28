classdef Rocket < handle

    properties
        name (1,1) string
        totalLength (1,1) double
        componentList dictionary
        aeroData

    end

    methods
        function obj = Rocket(name)
            obj.name = name;
            obj.componentList = dictionary();
        end


        function addComponent(rocketObj, componentObj)
            arguments
                rocketObj Rocket
                componentObj RocketComponent
            end

            numComponents = numEntries(rocketObj.componentList);

            if numComponents == 0
                rocketObj.componentList(componentObj.name) = {componentObj};
            elseif isKey(rocketObj.componentList, componentObj.name)
                %this should be removed for GUI implementation, only added here for testing
                warning("Replacing existing component '%s', do you want to proceed? (Y/N)", componentObj.name)
                response = input("", "s");
                if strcmp(response, "Y")
                    rocketObj.componentList(componentObj.name) = {componentObj};
                else
                    fprintf("\nStopping . . .\n")
                end
            else
                rocketObj.componentList(componentObj.name) = {componentObj};
            end
        end


        function removeComponent(rocketObj, componentName)
            if iskey(rocketObj.componentList, componentName)

                rocketObj.componentList = remove(rocketObj.componentList, componentName);
            else
                %this should be removed for GUI implementation, only added here for testing
                warning("The component you are trying to remove does not exist, try using a different name")
            end
        end


        function m = getMass(rocketObj)
            arguments
                rocketObj Rocket
            end

            componentArray = values(rocketObj.componentList);
            mass = cellfun(@(c) c.mass, componentArray);
            m = sum(mass);
        end


        function set.aeroData(rocketObj, filename)
            arguments
                rocketObj Rocket
                filename (1,1) string
            end

            filepath = "TheSixDoF" + filesep + "Inputs" + filesep + "RASAero" + filesep + filename + ".csv";
            rawData = readmatrix(filepath);
            data = [[rawData(1:300,1:5) rawData(1:300,8) rawData(1:300,13:15)]; 
                    [rawData(2501:2800,1:5) rawData(2501:2800,8) rawData(2501:2801,13:15)]; 
                    [rawData(5001:5300,1:5) rawData(5001:5300,8) rawData(5001:5301,13:15)]];
            rocketObj.aeroData = data;

        end

        function saveRocket(rocketObj)
            filename = rocketObj.name;
            filepath = "TheSixDoF" + filesep + "Inputs" + filesep + "Saved Rockets" + filesep + filename + ".mat";
            save(filepath, "rocketObj")
        end

    end
end