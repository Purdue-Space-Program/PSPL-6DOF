classdef Rocket < handle

    properties
        name (1,1) string
        componentArray cell
        componentDict dictionary
        %aeroData path

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
                rocketObj.componentArray{numComponents+1} = {componentObj};
                rocketObj.componentDict(componentObj.name) = componentObj;
            elseif isKey(rocketObj.componentDict, componentObj.name)
                warning("Replacing existing component '%s', do you want to proceed? (Y/N)", componentObj.name)
                response = input("", "s");
                if strcmp(response, "Y")
                    rocketObj.componentArray{numComponents+1} = {componentObj};
                    rocketObj.componentDict(componentObj.name) = componentObj;
                else
                    fprintf("\nStopping . . .\n")
                end
            else
                rocketObj.componentArray{numComponents+1} = {componentObj};
                rocketObj.componentDict(componentObj.name) = componentObj;
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

        function cm = getCoM(rocketObj)
            arguments
                rocketObj Rocket
            end

            cm = 0;
            m = rocketObj.getMass;

            for idx = length(rocketObj.componentArray)
                currentComponent = rocketObj.componentArray{idx};
                cm = cm + (currentComponent{1}.position * currentComponent{1}.mass);
            end

            cm = cm / m;
        end

    end
end