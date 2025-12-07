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
% Improve structural CoM update to include nosecone, struts, etc.
% Write MoI calcs
% Finish some smaller components (Avi, etc.) w/point mass inclusion

%---Initializations--------------------------------------------------------
compList = rocket.ComponentList;        % Write component list
numComponents = numEntries(compList);   % Total number of components
rocketMass = rocket.TotalMass;          % Total mass of the rocket [kg]
lengthTot = rocket.TotalLength;              % Total length of rocket [m]
radiusTot = rocket.OuterDiameter/2;     % Total radius of rocket [m]

%---CoM/MoI----------------------------------------------------------------
% Recover Propulsion Characteristics
if numComponents > 0
    values = compList.values;

    % Propulsion Properties
    for i = 1:numComponents
        if isa(values{i},'PropulsionSystem')
            propType = 'Liquid';
            propSysObj = values{i};
            burnTime = propSysObj.BurnTime;

        elseif isa(values{i},'SolidMotor')
            propType = 'Solid';
            propSysObj = values{i};
            burnTime = propSysObj.BurnTime;

        end
    end

    timeSpan = linspace(0,burnTime,burnTime*100)'; % Array for burn time
    numelTime = burnTime*100;     % Number of elements of time array
    CoMX = zeros(numelTime,1);
    CoMY = CoMX;
    CoMZ = CoMX;
    CoM = [CoMX, CoMY, CoMZ, timeSpan];    % [x|y|z|t] CoM position|time
    rocketMass = rocketMass .* ones(numelTime,1);
    countedMass = zeros(numelTime,1);      % Counted mass wrt time

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

        % Solid motor
        elseif isa(values{idx},'SolidMotor')
            
            % Object Values
            motorObj = values{i};
            burnTime = motorObj.BurnTime;
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
            CoMX = (CoM(:,1).*countedMass+CoMMotorX.*propellantMass)./(countedMass+propellantMass);
            CoMY = (CoM(:,2).*countedMass+CoMMotorY.*propellantMass)./(countedMass+propellantMass);
            CoMZ = (CoM(:,3).*countedMass+CoMMotorZ.*propellantMass)./(countedMass+propellantMass);
            CoM = [CoMX, CoMY, CoMZ];
            countedMass = countedMass + propellantMass;

        % Fins
        elseif isa(values{idx},'Fins')

            % Object values
            finObj = values{idx};
            rootChord = finObj.RootChord;
            tipChord = finObj.TipChord;
            mass = finObj.Mass;

            finMass = mass.*ones(numelTime,1);
            countedMass = countedMass + finMass;

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

        end
    end

    % Calculate Structural CoM/MoI
    % Note: currently assuming completely uniform mass distribution accross a
    % cylinder. Will change this in the future to improve accuracy by modeling
    % nosecone, individual struts, etc.

    % Structure CoM
    structureMass = rocketMass - countedMass(1);
    CoMStructX = (lengthTot/2).*(ones(numelTime,1));
    CoMStructY = 0;
    CoMStructZ = 0;

    % Total CoM Update
    CoMX = (CoM(:,1).*countedMass+CoMStructX.*structureMass)./(structureMass+countedMass);
    CoMY = (CoM(:,2).*countedMass+CoMStructY.*structureMass)./(structureMass+countedMass);
    CoMZ = (CoM(:,3).*countedMass+CoMStructZ.*structureMass)./(structureMass+countedMass);
    CoM = [CoMX, CoMY, CoMZ, timeSpan];

end

end