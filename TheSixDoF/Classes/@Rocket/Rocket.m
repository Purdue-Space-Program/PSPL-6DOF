classdef Rocket < handle

    properties
        name (1,1) string
        componentArray cell
        componentMap containers.Map

    end

    methods
        function obj = Rocket(name)
            obj.name = name;
            obj.componentArray = {};
            obj.componentMap = containers.Map("Keytype", 'char', "ValueType", 'any');
        end

        function addComponent(rocketObj, componentObj)
            if isKey(rocketObj.componentMap, componentObj.name)
                warning("Replacing existing component '%c', do you want to proceed? (Y/N)", componentObj.name)
                response = input("", "s");
                if response == ("Y" | "y" | "yes" | "Yes")
                    fprintf("\nStopping . . .\n")
                else
                    rocketObj.componentArray(end+1) = componentObj;
                    rocketObj.componentMap(componentObj.name) = componentObj;
                end
            else
                rocketObj.componentArray(end+1) = componentObj;
                rocketObj.componentMap(componentObj.name) = componentObj;
            end
        end

        function m = getMass(Rocket)
            m = 0;

            for idx = length(Rocket.componentArray)
                m = m + Rocket.componentArray(idx).mass;
            end
        end

        function cm = getCoM(Rocket)
            cm = 0;
            m = Rocket.getMass;

            for idx = length(Rocket.componentArray)
                cm = cm + (Rocket.componentArray(idx).position * Rocket.componentArray(idx).mass);
            end

            cm = cm / m;
        end

    end
end