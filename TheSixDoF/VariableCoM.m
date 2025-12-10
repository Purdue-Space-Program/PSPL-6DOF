function [totCoM, totMass, MoI] = VariableCoM(dt, tspan, graph, rocket)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSP FLIGHT DYNAMICS:
%
% Title: VariableCoM
% Author: Preston Wright - Created: 9/28/2024
%
% Description: This function calculates an array for the center of mass in 
% meters from the nose of the rocket. Uses given initialized parameters for 
% the rocket to approximate the center of mass 
%
% Inputs: 
% dt = given simulation time step [s]
% tspan = array of time values for total simulation run time for the given
%         time step [s]
% graph = boolean operator that controls the output of CoM visualizations:
%         1 will output visuals, 0 will not
%
% Outputs:
% totCoM  = array of (time, x coordinate of center of mass)
% totMass = array of (time, total mass)
% MoI     = array of (3x3 moment of inertia tensor)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Individual testing
    % Comment the function header out, uncomment this for testing
    % Uncomment function header, comment this out for normal use
% dt = 0.1;
% time = 13;
% arrayLength = (time / dt);
% tspan = linspace(0,time,arrayLength+1);
% graph = 1;

%% Conversions

% Initialize all the conversion multipliers
% NOTE: all conversions are scalar multiples in the code
IN2M = 0.0254;                          % Inches to meters
LBM2KG = 0.45359237;                    % Pound mass to kilograms
FT2M = 0.3048;                          % Feet to meters
IN2M3 = 0.00001639;                     % Inches cubed to meters cubed
densImp2Met = 27679.9;                  % Imperial to metric density units


%% Initializations

% Propellant initializations
massDry = 35.85782075197;                       % Initial dry mass [kg]
massWetInit = 39.94015208391423;                   % Initial wet mass [kg]
ullageFactor = 0.95;                    % Tank fill factor 
loxVol = 118.3861 * ullageFactor;        % Initial lox vol [in^3]
fuelVol = 155.6105 * ullageFactor;       % Initial fuel vol [in^3]
loxDens = 0.03922015;                   % Density of lox [kg/m^3]
fuelDens = 0.01450439;                  % Density of fuel [kg/m^3]

% Flow rate initializations
loxFlow = 1.91 * LBM2KG;               % Flow rate of lox [kg/s]
fuelFlow = 1.91 * LBM2KG;              % Flow rate of fuel [kg/s]
rocketHeight = 8.26 * FT2M;            % Total height of rocket [m]
tankResidual = 0.1;                    % Amount leftover 

% Tank size initializations
tankOD = 6 * IN2M;                  % Outer diameter of tank [m]
wallThick = 0.125 * IN2M;               % Wall thickness [m]
tankID = tankOD - 2 * wallThick;        % Inner diameter of tank [m]
tankIArea = pi * (tankID/2)^2;          % Inner tank area [m^2]

% Mass initializations
noseMass = 1.53;                 % Nose mass [kg]
heTMass = 6;                  % He tank mass [kg]
midAFMass = 0.5;                % Inner stage mass [kg]
empLoxTMass = 2.545;            % Empty lox tank mass [kg]
empFuelTMass = 2.695;           % Empty fuel tank mass [kg]
finCanMass = 9;               % Fin can mass [kg]
engineMass = 6;               % Engine mass [kg]

% Height initializations
noseHeight = 15 * IN2M;                 % Height of the nose [m]
heHeight = 18 * IN2M;                   % He tank height [m]
midAFHeight = 5 * IN2M;                % Inner stage height [m]
loxTHeight = 7 * IN2M;               % Lox tank height [m]
fuelTHeight = 8.44 * IN2M;                % Fuel tank height [m]
finCanHeight = 12 * IN2M;               % Fin can height [m]
engineHeight = 10.65 * IN2M;            % Engine height [m]
totHeight = noseHeight + heHeight ...
    + midAFHeight + loxTHeight + ...
    fuelTHeight + finCanHeight + ...
    engineHeight;                       % Total rocket height [m]

% Array initializations
loxRelCoM = zeros(length(tspan),1);     % Relative CoM for lox in tank [m]
loxMassArr = zeros(length(tspan),1);    % Lox total mass [m]
fuelRelCoM = zeros(length(tspan),1);    % Relative CoM for fuel in tank [m]
fuelMassArr = zeros(length(tspan),1);   % Fuel total mass [m]
loxCoM = zeros(length(tspan),1);        % Lox CoM [m]
fuelCoM = zeros(length(tspan),1);       % Fuel CoM [m]
totCoM = zeros(length(tspan),2);        % Total CoM [m]
totMass = zeros(length(tspan),2);       % Total Mass [m]
MoI = zeros(length(tspan),3,3);         % Moment of Inertia

%% Initial Calculations

% Initialized heights from nose
noseHFore = noseHeight/2;
heHFore = noseHeight + heHeight/2;
midAFHFore = noseHeight + heHeight + midAFHeight/2;
loxTHFore = noseHeight + heHeight + midAFHeight + loxTHeight/2;
fuelTHFore = noseHeight + heHeight + midAFHeight + loxTHeight + ...
    fuelTHeight/2;
finCanHFore = noseHeight + heHeight + midAFHeight + loxTHeight + ...
    fuelTHeight + finCanHeight/2;
engineHFore = noseHeight + heHeight + midAFHeight + loxTHeight + ...
    fuelTHeight + finCanHeight + engineHeight/2;

% Calculate empty CoM measured from nose
empCoM = (engineMass * engineHFore + finCanMass * finCanHFore + empFuelTMass ...
          * fuelTHFore + empLoxTMass * loxTHFore + midAFMass * midAFHFore ...
          + heTMass * heHFore + noseMass * noseHFore) / massDry;

% Measure height to fuel and lox from nose
heightToLox = noseHeight + heHeight + midAFHeight;
heightToFuel = heightToLox + loxTHeight;

% Measure propellant masses/CoM's from top of respective tank
loxMass = loxVol * loxDens * LBM2KG;
fuelMass = fuelVol * fuelDens * LBM2KG;
loxInitCoM = loxTHeight - (loxVol * ullageFactor * IN2M3 / tankIArea / 2);
fuelInitCoM = fuelTHeight - (fuelVol * ullageFactor * IN2M3 / tankIArea / 2);
finalLoxMass = loxMass / ullageFactor * tankResidual;
finalFuelMass = fuelMass / ullageFactor * tankResidual;

%% Calculations

% Loop through entire time array, updating CoM values along the way
for i = 1:length(tspan)

    % Initialize the first index in every array
    if i == 1

        % Start with initalized values above
        loxRelCoM(i) = loxInitCoM;
        loxMassArr(i) = loxMass;
        fuelRelCoM(i) = fuelInitCoM;
        fuelMassArr(i) = fuelMass;
        
        % Find CoM for propellants by adding relative CoM to tank location
        loxCoM(i) = loxRelCoM(i) + heightToLox;
        fuelCoM(i) = fuelRelCoM(i) + heightToFuel;

        % Initialize total CoM
        totCoM(i,2) = (loxCoM(i) * loxMassArr(i) + fuelCoM(i) * fuelMassArr(i) ...
                      + empCoM * massDry) / (loxMassArr(i) + fuelMassArr(i) ...
                      + massDry);
        totMass(i,2) = loxMassArr(i) + fuelMassArr(i) + massDry;
        
        % Update total simulation time
        totCoM(i,1) = tspan(i);
        totMass(i,1) = tspan(i);

    % Continue for all timesteps where both propellant masses are above the
    % final calculated mass (should approximately coincide with burn time)
    elseif loxMassArr(i-1)>finalLoxMass && fuelMassArr(i-1)>finalFuelMass
        
        % Update propellant masses
        loxMassArr(i) = loxMassArr(i-1) - loxFlow * dt;
        fuelMassArr(i) = fuelMassArr(i-1) - fuelFlow * dt;

        % Update relative propellant CoMs
        loxRelCoM(i) = loxTHeight - (loxMassArr(i) / (loxDens * ...
            densImp2Met) / tankIArea / 2);
        fuelRelCoM(i) = fuelTHeight - (fuelMassArr(i) / (fuelDens * ...
            densImp2Met) / tankIArea / 2);

        % Update propellant CoMs from nose
        loxCoM(i) = loxRelCoM(i) + heightToLox;
        fuelCoM(i) = fuelRelCoM(i) + heightToFuel;
        
        % Calculate total CoM from nose
        totCoM(i,2) = (loxCoM(i) * loxMassArr(i) + fuelCoM(i) * fuelMassArr(i) ...
                      + empCoM * massDry) / (loxMassArr(i) + fuelMassArr(i) ...
                      + massDry);
        totMass(i,2) = loxMassArr(i) + fuelMassArr(i) + massDry;
        
        % Update total simulation time
        totCoM(i,1) = tspan(i);
        totMass(i,1) = tspan(i);

    % Continue if fuel flow is zero (should approximately be post burn
    % time)
    else

        % Use final propellant masses
        loxMassArr(i) = finalLoxMass;
        fuelMassArr(i) = finalFuelMass;

        % Update relative propellant CoMs 
        loxRelCoM(i) = loxTHeight - (loxMassArr(i) / (loxDens * ...
            densImp2Met) / tankIArea / 2);
        fuelRelCoM(i) = fuelTHeight - (fuelMassArr(i) / (fuelDens * ...
            densImp2Met) / tankIArea / 2);

        % Update propellant CoMs from nose
        loxCoM(i) = loxRelCoM(i) + heightToLox;
        fuelCoM(i) = fuelRelCoM(i) + heightToFuel;
        
        % Update total CoM from nose
        totCoM(i,2) = (loxCoM(i) * loxMassArr(i) + fuelCoM(i) * fuelMassArr(i) ...
                      + empCoM * massDry) / (loxMassArr(i) + fuelMassArr(i) ...
                      + massDry);
        totMass(i,2) = loxMassArr(i) + fuelMassArr(i) + massDry;

        % Update total simulation time
        totCoM(i,1) = tspan(i);
        totMass(i,1) = tspan(i);
    end
    % Update inertia
    com = totCoM(i,2);
    MoI(i,1,1) = 1/2 * totMass(i,2) * (tankOD/2)^2;
    transverseMoI = 0;
    transverseMoI = transverseMoI + loxMassArr(i) * (1/4*(tankOD/2)^2 + 1/12*loxTHeight^2 + (loxCoM(i)-com)^2);
    transverseMoI = transverseMoI + fuelMassArr(i) * (1/4*(tankOD/2)^2 + 1/12*fuelTHeight^2 + (fuelCoM(i)-com)^2);
    transverseMoI = transverseMoI + engineMass * (1/4*(tankOD/2)^2 + 1/3*engineHeight^2 + (engineHFore-com)^2);
    transverseMoI = transverseMoI + finCanMass * (1/4*(tankOD/2)^2 + 1/3*finCanHeight^2 + (finCanHFore-com)^2);
    transverseMoI = transverseMoI + midAFMass * (1/4*(tankOD/2)^2 + 1/3*midAFHeight^2 + (midAFHFore-com)^2);
    transverseMoI = transverseMoI + heTMass * (1/4*(tankOD/2)^2 + 1/3*heHeight^2 + (heHFore-com)^2);
    %the nose is currently a cylinder
    transverseMoI = transverseMoI + noseMass * (1/4*(tankOD/2)^2 + 1/3*noseHeight^2 + (noseHFore-com)^2);
    
    MoI(i,2,2) = transverseMoI;
    MoI(i,3,3) = transverseMoI;
end

%% Optional CoM graphing output
if graph == 1

    % Initialize arrays and text used for graphing
    totCoMTop = totCoM(:,2);
    loxRelCoMIndividual = loxTHeight - loxRelCoM(:,1);
    fuelRelCoMIndividual = fuelTHeight - fuelRelCoM(:,1);
    txtRock = ['Rocket height: ' num2str(totHeight) ' meters'];
    txtLox = ['Lox tank height: ' num2str(loxTHeight) ' meters'];
    txtFuel = ['Fuel tank height: ' num2str(fuelTHeight) ' meters'];

    % Figure 1 - total CoM of the rocket
    figure(1)
    plot(tspan, totCoMTop, lineWidth=1.5)
    text(25,2.975,txtRock)
    grid on
    title("totalCoM versus time")
    xlabel("time")
    ylabel("CoM from nose [m]")

    % Figure 2 - CoM of the lox tank
    figure(2)
    plot(tspan, loxRelCoM, lineWidth=1.5)
    hold on
    plot(tspan, loxRelCoMIndividual, lineWidth=1.5)
    text(25,0.45,txtLox)
    grid on
    title("loxCoM from top of tank")
    xlabel("time")
    ylabel("loxCoM from top of lox tank")
    legend("From top", "From bottom", location="best")

    % Figure 3 - CoM of the fuel tank
    figure(3)
    plot(tspan, fuelRelCoM, lineWidth=1.5)
    hold on
    plot(tspan, fuelRelCoMIndividual, lineWidth=1.5)
   
    text(25,0.45,txtFuel)
    grid on
    title("fuelCoM from top of tank")
    xlabel("time")
    ylabel("fuelCoM from top of fuel tank")
    legend("From top", "From bottom", location="best")
    hold off

    figure(4)
    plot(tspan, MoI(:,1,1))
    hold on
    title('Dynamic MoI values')
    xlabel('time [s]')
    ylabel('MoI $I_{xx}$ [$kg \cdot m^2$]')

    yyaxis right
    
    plot(tspan,MoI(:,2,2))
    ylabel('MoI $I_{yy}$ [$kg \cdot m^2$]')


    
    legend('$I_{xx}$','$I_{yy}$/$I_{zz}$')
end
