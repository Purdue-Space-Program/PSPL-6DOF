function [totMass, CoM, MoI, MoIDot] = VariableCoMMoI(rocket)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSP FLIGHT DYNAMICS:
%
% Title: VariableCoMMoI
% Author: Preston Wright - Created: 11/9/2025
%
% Description: This takes in a rocket object with all necessary properties,
% and calculates the center of mass and inertia tensor with respect to
% time, from ignition through the end of burn.
%
% Inputs: 
% rocket = rocket object to determine the center of mass and inertia tensor
% for
%
% Outputs:
% totMass = Total mass of the vehicle wrt time [m|t] [kg|s]
% CoM = Center of mass of the vehicle as distance frome nose [x|y|z|t] [m]
% MoI = Inertia tensor of the vehicle [x|y|z|t] [m^4]
% MoIDot = Time derivative of inertia tensor [x|y|z|t] [m^4/s]
%
% Notes: All column vectors. Assume PMOI system. Assume CoM lies along
% lingitudinal axis of the rocket.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Future Updates:
% Improve structural CoM update to include nosecone, struts, etc.
% Validate MoI calcs
% Finish some smaller components (Avi, etc.) w/point mass inclusion

%---Initializations--------------------------------------------------------
compList = rocket.ComponentList;        % Write component list
numComponents = numEntries(compList);   % Total number of components
rocketMass = rocket.TotalMass;          % Total mass of the rocket [kg]
lengthTot = rocket.TotalLength;         % Total length of rocket [m]
radiusTot = rocket.OuterDiameter/2;     % Total radius of rocket [m]

%---CoM--------------------------------------------------------------------
% Calculations for rockets with components
if numComponents > 0
    values = compList.values;

    % Propulsion properties - burn time
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

    % Specific initializations
    timeSpan = linspace(0,burnTime,burnTime*100)'; % Array for burn time
    numelTime = burnTime*100;     % Number of elements of time array
    CoMX = zeros(numelTime,1);
    CoMY = CoMX;
    CoMZ = CoMX;
    MoIX = zeros(numelTime,1);
    MoIY = MoIX;
    MoIZ = MoIX;
    MoIXDot = zeros(numelTime,1);
    MoIYDot = MoIXDot;
    MoIZDot = MoIXDot;
    CoM = [CoMX, CoMY, CoMZ, timeSpan];    % [x|y|z|t] CoM position|time
    MoI = [MoIX, MoIY, MoIZ, timeSpan];    % [x|y|z|t] MoI value|time
    MoIDot = [MoIXDot, MoIYDot, MoIZDot, timeSpan];
    rocketMass = rocketMass .* ones(numelTime,1);
    countedMass = zeros(numelTime,1);      % Counted mass wrt time

    % Component-wise mass updates
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
            tankX = (position(1) + (length/2)).*ones(numelTime,1);
            tankY = position(2).*ones(numelTime,1);
            tankZ = position(3).*ones(numelTime,1);

            % Fluid System: assume constant mass flow, no slosh, no ullage
            % (Will likely change drain rate to use fuel/ox flow later with
            % ullage implementation)
            drainRateHeight = length/burnTime; 
            drainRateMass = massLiquid/burnTime;
            liquidLevel = length - drainRateHeight.*timeSpan;
            liquidMass = massLiquid - drainRateMass.*timeSpan;
            liquidLocate = length - liquidLevel;

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
            CoM = [CoMX, CoMY, CoMZ, timeSpan];
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
            CoM = [CoMX, CoMY, CoMZ, timeSpan];
            countedMass = countedMass + propellantMass;

            % MoI Update

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

    % Calculate Structural CoM
    % Note: currently assuming completely uniform mass distribution accross a
    % cylinder. Will change this in the future to improve accuracy by modeling
    % nosecone, individual struts, etc.

    % Structure CoM
    structureMass = rocketMass - countedMass(1);
    if structureMass < 0
        structureMass = 0;
    end
    CoMStructX = (lengthTot/2).*(ones(numelTime,1));
    CoMStructY = 0;
    CoMStructZ = 0;

    % Final total CoM update
    CoMX = (CoM(:,1).*countedMass+CoMStructX.*structureMass)./(structureMass+countedMass);
    CoMY = (CoM(:,2).*countedMass+CoMStructY.*structureMass)./(structureMass+countedMass);
    CoMZ = (CoM(:,3).*countedMass+CoMStructZ.*structureMass)./(structureMass+countedMass);
    CoM = [CoMX, CoMY, CoMZ, timeSpan];

    totMass = [structureMass + countedMass, timeSpan];

%---MoI--------------------------------------------------------------------

    % Component-wise inertia updates
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

            % Tank System
            tankX = (position(1) + (length/2)).*ones(numelTime,1);
            tankY = position(2).*ones(numelTime,1);
            tankZ = position(3).*ones(numelTime,1);

            % Fluid System: assume constant mass flow, no slosh, no ullage
            % (Will likely change drain rate to use fuel/ox flow later with
            % ullage implementation)
            drainRateHeight = length/burnTime; 
            drainRateMass = massLiquid/burnTime;
            liquidLevel = length - drainRateHeight.*timeSpan;
            liquidMass = massLiquid - drainRateMass.*timeSpan;
            liquidLocate = length - liquidLevel;
            
            liquidX = position(1)+liquidLocate+(liquidLevel./2);
            liquidY = tankY;
            liquidZ = tankZ;

            % Dry tank inertia
            tankMoIX = 0.5*mass*(radius^2+(radius-thick)^2);
            tankMoIY = (1/12)*mass*(3*(radius^2+(radius-thick)^2)+length^2);
            tankMoIZ = tankMoIY;
                                 
            % Dry tank parallel axis
            MoIX = MoI(:,1) + (tankMoIX + mass.*(CoM(:,2).^2+CoM(:,3).^2));
            MoIY = MoI(:,2) + (tankMoIY + mass.*(abs(tankX-CoM(:,1)).^2+CoM(:,3).^2));
            MoIZ = MoI(:,3) + (tankMoIZ + mass.*(abs(tankX-CoM(:,1)).^2+CoM(:,2).^2));
            MoI = [MoIX, MoIY, MoIZ, timeSpan];

            % Liquid inertia
            LiquidMoIX = 0.5.*liquidMass.*(radius-thick).^2;
            LiquidMoIY = (1/12).*liquidMass.*(3.*(radius-thick).^2+liquidLevel.^2);
            LiquidMoIZ = LiquidMoIY;

            % Liquid parallel axis
            MoIX = MoI(:,1) + (LiquidMoIX + liquidMass.*(CoM(:,2).^2+CoM(:,3).^2));
            MoIY = MoI(:,2) + (LiquidMoIY + liquidMass.*(abs(liquidX-CoM(:,1)).^2+CoM(:,3).^2));
            MoIZ = MoI(:,3) + (LiquidMoIZ + liquidMass.*(abs(liquidX-CoM(:,1)).^2+CoM(:,2).^2));
            MoI = [MoIX, MoIY, MoIZ, timeSpan];

        % Solid motor
        elseif isa(values{idx},'SolidMotor')
            
            % Object Values
            motorObj = values{i};
            burnTime = motorObj.BurnTime;
            mass = motorObj.Mass;
            position = motorObj.Position;
            length = motorObj.Length;

            motorX = position(1) + propellantLength/2;
            motorY = position(2).*ones(numelTime,1);
            motorZ = position(3).*ones(numelTime,1);

            % Propellant System
            % (Will change later to better model solid motor burning (input
            % propellant geometry))
            engineLengthRate = length/burnTime;
            massLossRate = mass/burnTime;
            propellantLength = length - engineLengthRate.*timeSpan;
            propellantMass = mass - massLossRate.*timeSpan;
            propellantLocate = length - propellantLength;
            propellantX = position(1)+propellantLocate+(propellantLength./2);

            % Solid motor inertia
            MotorMoIX = 0.5.*propellantMass.*(radiusTot^2);
            MotorMoIY = (1/12).*propellantMass.*(3*(radiusTot)^2+propellantLength.^2);
            MotorMoIZ = MotorMoIY;

            % Solid motor parallel axis
            MoIX = MoI(:,1) + (MotorMoIX + propellantMass.*(CoM(:,2).^2+CoM(:,3).^2));
            MoIY = MoI(:,2) + (MotorMoIY + propellantMass.*(abs(CoM(:,1)-propellantX).^2+CoM(:,3).^2));
            MoIZ = MoI(:,3) + (MotorMoIZ + propellantMass.*(abs(CoM(:,1)-propellantX).^2+CoM(:,2).^2));
            MoI = [MoIX, MoIY, MoIZ, timeSpan];

        end
    end

% Inertia calculation for structure
structMoIX = 0.5*structureMass*radiusTot^2;
structMoIY = (1/12)*structureMass*(3*radiusTot^2+lengthTot^2);
structMoIZ = structMoIY;

% Parallel axis for structure/final MoI update
MoIX = MoI(:,1) + (structMoIX + structureMass.*(CoM(:,2).^2+CoM(:,3).^2));
MoIY = MoI(:,2) + (structMoIY + structureMass.*(abs(lengthTot/2-CoM(:,1)).^2+CoM(:,3).^2));
MoIZ = MoI(:,3) + (structMoIZ + structureMass.*(abs(lengthTot/2-CoM(:,1)).^2+CoM(:,2).^2));
MoI = [MoIX, MoIY, MoIZ, timeSpan];

% Calculations for rocket without components
elseif numComponents == 0
    
    % Use CoM, MoI of a cylinder with uniform density
    % CoM Calcs
    CoMX = lengthTot/2;
    CoMY = 0;
    CoMZ = 0;
    CoM = [CoMX, CoMY, CoMZ];

    % MoI Calcs
    MoIX = (1/2)*rocketMass*radiusTot^2;
    MoIY = (1/12)*rocketMass*(3*radiusTot^2+lengthTot^2);
    MoIZ = (1/12)*rocketMass*(3*radiusTot^2+lengthTot^2);
    MoI = [MoIX, MoIY, MoIZ];

end
end