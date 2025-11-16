function [CoM, MoI] = VariableCoMMoI(rocket)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSP FLIGHT DYNAMICS:
%
% Title: VariableCoMMoI
% Author: Preston Wright - Created: 11/9/2025
%
% Description: This takes in a rocket object with all necessary properties,
% and calculates the center of mass and inertia tensor with respect to
% time, from ignition through the end of burn. Note: all column vectors.
%
% Inputs: 
% rocket = rocket object to determine the center of mass and inertia tensor
% for
%
% Outputs:
% CoM = Center of mass of the vehicle as distance frome nose [x|y|z|t] [m]
% MoI = Inertia tensor of the vehicle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Future Updates:
% Actually test (my CMS file isn't downloading properly)
% Include structural CoM (need overall object update)
% Write MoI calcs
% Finish some smaller components

%---Initializations--------------------------------------------------------
compList = rocket.ComponentList;        % Write component list
numComponents = numEntries(compList);   % Total number of components
lengthTot = rocket.Length;              % Total length of rocket [m]
radiusTot = rocket.OuterDiameter/2;     % Total radius of rocket [m]

%---Individual CoM/MoI-----------------------------------------------------
% Recover individual CoM/MoI values of components
if numComponents > 0
    values = compList.values;

    % Propulsion Properties
    for i = 1:numComponents
        if isa(values{i},'PropulsionSystem')
            propType = 'Liquid';
            propSysObj = values{i};
            burnTime = propSysObj.BurnTime;
            fuelFlow = propSysObj.FuelMassFlow;
            oxFlow = propSysObj.OxMassFlow;

        elseif isa(values{i},'SolidMotor')
            propType = 'Solid';
            propSysObj = values{i};
            burnTime = propSysObj.BurnTime;
            massFlow = propSysObj.MassFlow;

        end
    end

    timeSpan = linspace(0,burnTime,burnTime*100)'; % Array for burn time
    numelTime = length(timeSpan);     % Number of elements of time array
    CoMX = zeros(numelTime,1);
    CoMY = CoMX;
    CoMZ = CoMX;
    CoM = [CoMX, CoMY, CoMZ, timeSpan];    % [x|y|z|t] CoM position|time
    countedMass = zeros(numelTime,1);      % Counted mass wrt time

    % Will likely want precalc for just structural CoM

    % Component-wise mass/inertia updates
    for idx = 1:numComponents
        
        % Tanks
        if isa(values{idx},'Tank')
            
            % Object values
            tankObj = values{idx};
            mass = tankObj.Mass;
            massLiquid = tankObj.LiquidMass;
            length = tankObj.Length;
            radius = tankObj.TankDia/2;
            thick = tankObj.Thickness;
            contents = tankObj.FuelOx;
            position = tankObj.Position;

            % Tank CoM
            tankX = (position(1)+(length/2)).*ones(numelTime,1);
            tankY = position(2).*ones(numelTime,1);
            tankZ = position(3).*ones(numelTime,1);

            % Fluid System: assume constant mass flow, no slosh, no ullage
            % (Will likely change drain rate to use fuel/ox flow later with
            % ullage implementation)
            drainRateHeight = length/burnTime; 
            drainRateMass = massLiquid/burnTime;
            liquidLevel = (length-drainRateHeight.*timeSpan);
            liquidMass = massLiquid - drainRateMass.*timeSpan;
            liquidLocate = length-liquidLevel;

            % Liquid CoM
            liquidX = position(1)+liquidLocate+(liquidLevel/2);
            liquidY = tankY;
            liquidZ = tankZ;

            % Full Tank CoM
            tankTotX = (tankX.*mass+liquidX.*liquidMass)./(liquidMass+mass);
            tankTotY = (tankY.*mass+liquidY.*liquidMass)./(liquidMass+mass);
            tankTotZ = (tankZ.*mass+liquidZ.*liquidMass)./(liquidMass+mass);

            % Update CoM totals
            CoMX = (CoM(:,1).*countedMass+tankTotX.*(liquidMass+mass))./(countedMass+mass+liquidMass);
            CoMY = (CoM(:,2).*countedMass+tankTotY.*(liquidMass+mass))./(countedMass+mass+liquidMass);
            CoMZ = (CoM(:,3).*countedMass+tankTotZ.*(liquidMass+mass))./(countedMass+mass+liquidMass);
            CoM = [CoMX, CoMY, CoMZ];
            countedMass = countedMass + liquidMass + mass;

        elseif isa(values{idx},'SolidMotor')
            
            % Object Values
            motorObj = values{i};
            burnTime = motorObj.BurnTime;
            massFlow = motorObj.MassFlow;
            mass = motorObj.Mass;
            position = motorObj.Position;
            length = motorObj.Length;

            % Propellant System
            % (Will change later to better model solid motor burning (input
            % propellant geometry))
            engineLengthRate = length/burnTime;
            massLossRate = mass/burnTime;
            propellantLength = length - engineLengthRate.*timeSpan;
            propellantMass = mass - massLossRate.*timeSpan;

            % Motor CoM
            CoMMotorX = position(1) + propellantLength/2;
            CoMMotorY = position(2).*ones(numelTime,1);
            CoMMotorZ = position(3).*ones(numelTime,1);

            % Total CoM Update
            CoMX = (CoM(:,1).*countedMass+CoMMotorX.*propellantMass)/(countedMass+propellantMass);
            CoMY = (CoM(:,2).*countedMass+CoMMotorY.*propellantMass)/(countedMass+propellantMass);
            CoMZ = (CoM(:,3).*countedMass+CoMMotorZ.*propellantMass)/(countedMass+propellantMass);
            CoM = [CoMX, CoMY, CoMZ];
            countedMass = countedMass + propellantMass;

        % Fins
        elseif isa(values{idx},'Fins')

            % Object values
            finObj = values{idx};
            rootChord = finObj.RootChord;
            tipChord = finObj.TipChord;
            mass = finObj.Mass;

            % 


        %         Count (1,1) int8 {mustBeMember(Count, [3, 4])} = 3 % Fin Count
        % 
        % % add the properties for the airfoil later with mustBeMeber
        % 
        % Airfoil (1,1) string % Airfoil
        % Span (1,1) double % Span [m]
        % RootChord (1,1) double % Root Chord Length [m]
        % TipChord (1,1) double % Tip Chord Length [m]
        % Sweep (1,1) double % Sweep [m]
        % Thickness (1,1) double % Thickness [m]
         
        % Solid Motor

        end
    end
end