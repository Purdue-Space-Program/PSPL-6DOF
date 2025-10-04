classdef Rocket < handle

    properties
        name (1,1) string
        componentArray cell
        componentDict dictionary
        aeroData path

    end

    methods
        function obj = Rocket(name)
            obj.name = name;
            obj.componentArray = {};
            obj.componentDict = dictionary();
        end

        function addComponent(rocketObj, componentObj)
            arguments
                rocketObj Rocket
                componentObj RocketComponent
            end

            if isKey(rocketObj.componentDict, componentObj.name)
                warning("Replacing existing component '%c', do you want to proceed? (Y/N)", componentObj.name)
                response = input("", "s");
                if response == ("Y" | "y" | "yes" | "Yes")
                    fprintf("\nStopping . . .\n")
                else
                    rocketObj.componentArray(end+1) = componentObj;
                    rocketObj.componentDict(componentObj.name) = componentObj;
                end
            else
                rocketObj.componentArray(end+1) = componentObj;
                rocketObj.componentDict(componentObj.name) = componentObj;
            end
        end

        function m = getMass(rocketObj)
            arguments
                rocketObj Rocket
            end

            m = 0;

            for idx = length(rocketObj.componentArray)
                m = m + rocketObj.componentArray(idx).mass;
            end
        end

        function cm = getCoM(rocketObj)
            arguments
                rocketObj Rocket
            end

            cm = 0;
            m = rocketObj.getMass;

            for idx = length(rocketObj.componentArray)
                cm = cm + (rocketObj.componentArray(idx).position * rocketObj.componentArray(idx).mass);
            end

            cm = cm / m;
        end

    end
end