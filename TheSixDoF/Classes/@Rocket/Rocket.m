classdef Rocket < handle

    properties
        name (1,1) string
        totalLength (1,1) double
        componentArray cell
        componentDict dictionary
        aeroData

    end

    methods
        function obj = Rocket(name)
            obj.name = name;
            obj.componentArray = cell(1);
            obj.componentDict = dictionary();
        end

        function addComponent(rocketObj, componentObj)
            arguments
                rocketObj Rocket
                componentObj RocketComponent
            end

            numComponents = numEntries(rocketObj.componentDict);

            if numComponents == 0
                rocketObj.componentArray{numComponents+1} = componentObj;
                rocketObj.componentDict(componentObj.name) = componentObj;
            elseif isKey(rocketObj.componentDict, componentObj.name)
                %this should be removed for GUI implementation, only added here for testing
                warning("Replacing existing component '%s', do you want to proceed? (Y/N)", componentObj.name)
                response = input("", "s");
                if strcmp(response, "Y")
                    rocketObj.componentArray{numComponents+1} = componentObj;
                    rocketObj.componentDict(componentObj.name) = componentObj;
                else
                    fprintf("\nStopping . . .\n")
                end
            else
                rocketObj.componentArray{numComponents+1} = componentObj;
                rocketObj.componentDict(componentObj.name) = componentObj;
            end
        end

        function removeComponent(rocketObj, componentName)
            if iskey(rocketObj.componentDict, componentName)

                keyArray = keys(rocketObj.componentDict);
                componentIndex = strcmp(keyArray, componentName);

                rocketObj.componentArray(componentIndex) = [];
                rocketObj.componentDict = remove(rocketObj.componentDict, componentName);
            else
                %this should be removed for GUI implementation, only added here for testing
                warning("The component you are trying to remove does not exist, try using a different name")
            end
        end

        function m = getMass(rocketObj)
            arguments
                rocketObj Rocket
            end

            m = 0;

            for idx = length(rocketObj.componentArray)
                currentComponent = rocketObj.componentArray{idx};
                m = m + currentComponent{1}.mass;
            end
        end

        function set.aeroData(rocketObj, filename)
            arguments
                rocketObj Rocket
                filename (1,1) string
            end

            filepath = "\TheSixDoF\Inputs\RASAero\" + filename + ".csv";
            rawData = readmatrix(filepath);
            data = [[rawData(1:300,1:5) rawData(1:300,8) rawData(1:300,13:15)]; 
                    [rawData(2501:2800,1:5) rawData(2501:2800,8) rawData(2501:2801,13:15)]; 
                    [rawData(5001:5300,1:5) rawData(5001:5300,8) rawData(5001:5301,13:15)]];
            rocketObj.aeroData = data;

        end

    end
end