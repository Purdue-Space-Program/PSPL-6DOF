classdef Rocket
    % NOTE: the way the class constructor gets and stores rocket data will NOT work unless your root directory (in MATLAB or VSCode) is set to "PSPL-6DoF," so always check that before running

    properties
        name (1,1) string
        thrust (1,1) double {mustBePositive}
        burnTime (1,1) double {mustBePositive}
        exitArea (1,1) double {mustBePositive}
        exitPressure (1,1) double {mustBePositive}
        bodyDiameter (1,1) double {mustBePositive}
        totalLength (1,1) doule {mustBePositive}
        finCount (1,1) int {mustBePositive}
        finHeight (1,1) double {mustBePositive}
        finThickness (1,1) double {mustBePositive}
        finType (1,1) {mustBeMember(finType,['Delta', 'Swept', 'Trapezoidal', 'Irregular'])} = 'Trapezoidal'
        noseconeLength (1,1) double {mustBePositive}
        noseconeType (1,1) {mustBeMember(noseconeType,['Von Karman', 'Tangent Ogive', 'Conical', 'Parabolic', 'LV-Haack', 'LD-Haack', 'Power Series'])} = 'Power Series'
        wetMass (1,1) double {mustBePositive}
        dryMass (1,1) double {mustBePositive}
        fuel (1,1) string {mustBeMember(fuel,['LNG', 'Ethanol'])} = 'Ethanol'
        fuelMass (1,1) double {mustBePositive}
        fuelUllage (1,1) double {mustBePositive}
        fuelMassFlow (1,1) double {mustBePositive}
        fuelVolume (1,1) double {mustBePositive}
        fuelTankOuterDiameter (1,1) double {mustBePositive}
        fuelTankThickness (1,1) double {mustBePositive}
        oxidizer (1,1) string {mustBeMember(oxidizer, 'LOX')} = 'LOX'
        oxMass (1,1) double {mustBePositive}
        oxUllage (1,1) double {mustBePositive}
        oxMassFlow (1,1) double {mustBePositive}
        oxVolume (1,1) double {mustBePositive}
        oxTankOuterDiameter (1,1) double {mustBePositive}
        oxTankThickness (1,1) double {mustBePositive}
        injectorType (1,1) string {mustBeMember(injectorType,['Pintle', 'Impinging', 'Swirl', 'Shower Head'])} = 'Pintle'
        motorCoolingType (1,1) string {mustBeMember(motorCoolingType,['Heatsink', 'Regen', 'Ablative'])} = 'Ablative'
        tankPressMethod (1,1) string {mustBeMember(tankPressMethod,['Pressure-Fed', 'Pump-Fed'])} = 'Pressure-Fed'
        chamberPressure (1,1) double {mustBePositive}
        specificImpulse (1,1) double {mustBePositive}
        motorLength (1,1) double {mustBePositive}
        lowerAirframeLength (1,1) double {mustBePositive}
        fuelTankLength (1,1) double {mustBePositive}
        oxTankLength (1,1) double {mustBePositive}
        interTankLength (1,1) double {mustBePositive}
        upperAirframeLength (1,1) double {mustBePositive}
    end

    properties (Access = protected, Dependent)
        finWettedArea (1,1) double {mustBePositive}
        finFrontalArea (1,1) double {mustBePositive}
        aeroData (1,1) string
        finenessRatio (1,1) int {mustBeInRange(finenessRatio,1:10)}
        refArea (1,1) double {mustBePositive}
        fuelOxRatio (1,1) double {mustBePositive}
        copvMass (1,1) double {mustBePositive}
        copvVolume (1,1) double {mustBePositive}
        initialCoM (1,1) double {mustBePositive}
        finalCoM (1,1) double {mustBePositive}
    end

    methods

        function rocket = Rocket()
            name_prompt = "Input the name of a saved rocket (listed below), 1 to save a new rocket, or 2 to create a temporary object (not saved) -> ";
            name_struct = dir("TheSixDoF\Inputs\Saved Saved Rockets\*.mat");
            name_str = "";
            for index = 1:length(name_struct)
                current_name = name_struct(index).name;
                name_str = name_str + newline + convertCharsToStrings(current_name);
            end
            
            response = input(name_prompt + name_str,"s");
            if isempty(str2num(response))
                inputName = response;
                fprintf("\nSearching Saved Rockets . . .")
                nameValid = 0;
                while nameValid == 0
                    if find(inputName, name_str)
                        nameValid = 1;
                    else
                        fprintf("\nERROR: Rocket could not be found, please input a different name ->")
                    end
                end
                fprintf("\nRocket %s retrieved, loading saved data.", inputName)
                rocketDataPath = "TheSixDoF\Inputs\Saved Rockets\" + intputName + ".mat";

                rocketData = load(rocketDataPath);

                % Set all properties
                rocket.name = rocketData.name;

            elseif response == 1
                nameValid = 0;
                while nameValid == 0
                name_prompt = newline + "Please enter the name of the new rocket (note, the name can be changed later if needed) -> ";
                inputName = input(name_prompt, "s");
                    if find(inputName, name_str)
                        fprintf("\nERROR: Name is already in use, please select a new one.")
                    else
                        nameValid = 1;
                    end
                end

                rocketDataPath = "TheSixDoF\Inputs\Saved Rockets\" + inputName + ".mat";
                rocketData = matfile(rocketDataPath, "Writable", true);

                % Set all properties
                rocket.name = inputName;

                % Save all properties
                rocketData.name = inputName;

            elseif response == 2
                name_prompt = newline + "Please enter the name of the new rocket (note, the name can be changed later if needed) -> ";
                inputName = input(name_prompt, "s");

                % Set all properties
                rocket.name = inputName;
            end
        end

        function drawRocket(rocket)
            [xBody,yBody,zBody] = cylinder(rocket.bodyDiameter/2,50);

            figure()
            zBody = zBody * rocket.totalLength;
            surf(xBody,yBody,zBody, 'FaceColor','red', 'LineStyle','none', 'FaceAlpha','1')
            axis equal
        end

        % method for determining the area of a single fin parallel to the airstream
        function finWettedArea = get.finWettedArea(roc)
            if roc.finType == "Delta"
                prompt = newline + "Input fin wetted area in m^2 (Delta) -> ";
                finWettedArea = input(prompt);
            elseif roc.finType == "Trapezoidal"
                prompt = newline + "Input fin wetted area in m^2 (Trapezoidal) -> ";
                finWettedArea = input(prompt);
            elseif roc.finType == "Swept"
                prompt = newline + "Input fin wetted area in m^2 (Swept) -> ";
                finWettedArea = input(prompt);
            else
                prompt = newline + "Input fin wetted area in m^2 (Other) -> ";
                finWettedArea = input(prompt);
            end
        end

        % method for determining the area of a single fin normal to the airstream
        function finFrontalArea = get.finFrontalArea(roc)
            finFrontalArea = roc.finHeight * roc.finThickness;
        end

        % method for getting all RASAero data for the rocket
        % ideally making it dependent means we won't have to store it every time a rocket object is created,
        % but it might also slow things down, so this behavior could be changed later
        function aeroData = get.aeroData(roc)
            aeroDataPath = "TheSixDoF\Inputs\RASAero" + roc.name + ".csv";
            aeroData = readmatrix(aeroDataPath);
        end

        function finenessRatio = get.finenessRatio(roc)
            finenessRatio = roc.totalLength / roc.noseconeLength;
        end

        function refArea = get.refArea(roc)
            refArea = (pi * (roc.bodyDiameter / 2)^2) + (roc.finCount * roc.finFrontalArea);
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