function [totCoM, totMass, MoI] = VariableCoM(dt, tspan, plot_CoM_graph, rocket, rocket_name)
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

switch rocket_name
    case "Rocket A"
        
        % General
        rocket_diameter = 6 * IN2M; % OD of the rocket [in]
        rocket_height = 8.26 * FT2M;            % Total height of rocket [ft]


        % Propellant initializations
        mass_dry = rocket.TotalMass;             % Initial dry mass [kg]
        mass_wet_init = 27.93;        % Initial wet mass [kg]
        ullage_factor = 0.95;                    % Tank fill factor 
        lox_vol = 118.3861 * ullage_factor;       % Initial lox vol [in^3]
        fuel_vol = 155.6105 * ullage_factor;      % Initial fuel vol [in^3]
        LOx_density = 0.03922015;                   % Density of lox [kg/m^3]
        fuel_density = 0.01450439;                  % Density of fuel [kg/m^3]
        
        % Flow rate initializations
        lox_mass_flow_rate = 1.91 * LBM2KG;               % Flow rate of lox [kg/s]
        fuel_mass_flow_rate = 1.91 * LBM2KG;              % Flow rate of fuel [kg/s]
        tank_residual = 0.1;                    % Amount leftover 
        
        % Tank size initializations
        tank_OD = 6 * IN2M;                  % Outer diameter of tank [m]
        tank_wall_thickness = 0.125 * IN2M;           % Wall thickness [m]
        tank_ID = tank_OD - 2 * tank_wall_thickness;    % Inner diameter of tank [m]
        tank_inner_area = pi * (tank_ID/2)^2;      % Inner tank area [m^2]
        
        % Mass initializations
        nose_mass = 1.53;                 % Nose mass [kg]
        pressurant_tank_mass = 6;        % Pressurant tank mass [kg]
        mid_mass = 0.5;                 % Inner stage mass [kg]
        empty_LOx_tank_mass = 2.545;             % Empty lox tank mass [kg]
        empty_fuel_tank_mass = 2.695;            % Empty fuel tank mass [kg]
        fin_can_mass = 9;                  % Fin can mass [kg]
        engine_mass = 6;                  % Engine mass [kg]
        
        % Height initializations
        nose_height = 15 * IN2M;                 % Height of the nose [m]
        pressurant_tank_height = 18 * IN2M;     % Pressurant tank height [m]
        mid_height = 5 * IN2M;                 % Inner stage height [m]
        lox_tank_height = 7 * IN2M;                  % Lox tank height [m]
        fuel_tank_height = 8.44 * IN2M;              % Fuel tank height [m]
        fin_can_height = 12 * IN2M;               % Fin can height [m]
        engine_height = 10.65 * IN2M;            % Engine height [m]
    
        total_height = rocket_height;
    
    case "CMS"
        
        % Propellant initializations
        mass_dry = 49.877;                       % Initial dry mass [kg]
        mass_wet_init = 9999999999999;                   % Initial wet mass [kg]
        ullage_factor = 0.95;                    % Tank fill factor 
        lox_vol = 1201.10 * ullage_factor;        % Initial lox vol [in^3]
        fuel_vol = 1299.95 * ullage_factor;       % Initial fuel vol [in^3]
        LOx_density = 0.03922015;                   % Density of lox [kg/m^3]
        fuel_density = 0.01450439;                  % Density of fuel [kg/m^3]
        
        % Flow rate initializations
        lox_mass_flow_rate = 2.916 * LBM2KG;               % Flow rate of lox [kg/s]
        fuel_mass_flow_rate = 1.188 * LBM2KG;              % Flow rate of fuel [kg/s]
        rocket_height = 16.81 * FT2M;            % Total height of rocket [m]
        tank_residual = 0.07;                    % Amount leftover 
        
        % Tank size initializations
        tank_OD = 6.625 * IN2M;                  % Outer diameter of tank [m]
        tank_wall_thickness = 0.134 * IN2M;               % Wall thickness [m]
        tank_ID = tank_OD - 2 * tank_wall_thickness;        % Inner diameter of tank [m]
        tank_inner_area = pi * (tank_ID/2)^2;          % Inner tank area [m^2]
        
        % Mass initializations
        nose_mass = 11 * LBM2KG;                 % Nose mass [kg]
        pressurant_tank_mass = 18 * LBM2KG;     % Pressurant tank mass [kg]
        mid_mass = 23 * LBM2KG;                % Inner stage mass [kg]
        empty_LOx_tank_mass = 9.05 * LBM2KG;            % Empty lox tank mass [kg]
        empty_fuel_tank_mass = 9.91 * LBM2KG;           % Empty fuel tank mass [kg]
        fin_can_mass = 28 * LBM2KG;               % Fin can mass [kg]
        engine_mass = 11 * LBM2KG;               % Engine mass [kg]
        
        % Height initializations
        nose_height = 33 * IN2M;                 % Height of the nose [m]
        pressurant_tank_height = 43 * IN2M;     % Pressurant tank height [m]
        mid_height = 16 * IN2M;                % Inner stage height [m]
        lox_tank_height = 35.6 * IN2M;               % Lox tank height [m]
        fuel_tank_height = 39 * IN2M;                % Fuel tank height [m]
        fin_can_height = 22 * IN2M;               % Fin can height [m]
        engine_height = 11.77 * IN2M;            % Engine height [m]
        
        total_height = nose_height + pressurant_tank_height ...
            + mid_height + lox_tank_height + ...
            fuel_tank_height + fin_can_height + ...
            engine_height;                       % Total rocket height [m]
    
    otherwise
        error("Invalid rocket name!")
end



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
noseHFore = nose_height/2;
heHFore = nose_height + pressurant_tank_height/2;
midAFHFore = nose_height + pressurant_tank_height + mid_height/2;
loxTHFore = nose_height + pressurant_tank_height + mid_height + lox_tank_height/2;
fuelTHFore = nose_height + pressurant_tank_height + mid_height + lox_tank_height + ...
    fuel_tank_height/2;
finCanHFore = nose_height + pressurant_tank_height + mid_height + lox_tank_height + ...
    fuel_tank_height + fin_can_height/2;
engineHFore = nose_height + pressurant_tank_height + mid_height + lox_tank_height + ...
    fuel_tank_height + fin_can_height + engine_height/2;
        

if isprop(rocket, "CoMOverride") && ~isempty(rocket.CoMOverride)
    rocket_empty_CoM = rocket.CoMOverride; % [m]
    disp("exists!")
else
    % Calculate empty CoM measured from nose
    rocket_empty_CoM = (engineMass * engineHFore + finCanMass * finCanHFore + empFuelTMass ...
      * fuelTHFore + empLoxTMass * loxTHFore + midAFMass * midAFHFore ...
      + pressurant_tank_mass * heHFore + noseMass * noseHFore) / massDry;
    disp("doesnt exist")
end
fprintf("rocket_empty_CoM %f\n", rocket_empty_CoM)

% Measure height to fuel and lox from nose
heightToLox = nose_height + pressurant_tank_height + mid_height;
heightToFuel = heightToLox + lox_tank_height;

% Measure propellant masses/CoM's from top of respective tank
loxMass = lox_vol * LOx_density * LBM2KG;
fuelMass = fuel_vol * fuel_density * LBM2KG;
loxInitCoM = lox_tank_height - (lox_vol * ullage_factor * IN2M3 / tank_inner_area / 2);
fuelInitCoM = fuel_tank_height - (fuel_vol * ullage_factor * IN2M3 / tank_inner_area / 2);
finalLoxMass = loxMass / ullage_factor * tank_residual;
finalFuelMass = fuelMass / ullage_factor * tank_residual;

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
                      + rocket_empty_CoM * mass_dry) / (loxMassArr(i) + fuelMassArr(i) ...
                      + mass_dry);
        totMass(i,2) = loxMassArr(i) + fuelMassArr(i) + mass_dry;
        
        % Update total simulation time
        totCoM(i,1) = tspan(i);
        totMass(i,1) = tspan(i);

    % Continue for all timesteps where both propellant masses are above the
    % final calculated mass (should approximately coincide with burn time)
    elseif loxMassArr(i-1)>finalLoxMass && fuelMassArr(i-1)>finalFuelMass
        
        % Update propellant masses
        loxMassArr(i) = loxMassArr(i-1) - lox_mass_flow_rate * dt;
        fuelMassArr(i) = fuelMassArr(i-1) - fuel_mass_flow_rate * dt;

        % Update relative propellant CoMs
        loxRelCoM(i) = lox_tank_height - (loxMassArr(i) / (LOx_density * ...
            densImp2Met) / tank_inner_area / 2);
        fuelRelCoM(i) = fuel_tank_height - (fuelMassArr(i) / (fuel_density * ...
            densImp2Met) / tank_inner_area / 2);

        % Update propellant CoMs from nose
        loxCoM(i) = loxRelCoM(i) + heightToLox;
        fuelCoM(i) = fuelRelCoM(i) + heightToFuel;
        
        % Calculate total CoM from nose
        totCoM(i,2) = (loxCoM(i) * loxMassArr(i) + fuelCoM(i) * fuelMassArr(i) ...
                      + rocket_empty_CoM * mass_dry) / (loxMassArr(i) + fuelMassArr(i) ...
                      + mass_dry);
        totMass(i,2) = loxMassArr(i) + fuelMassArr(i) + mass_dry;
        
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
        loxRelCoM(i) = lox_tank_height - (loxMassArr(i) / (LOx_density * ...
            densImp2Met) / tank_inner_area / 2);
        fuelRelCoM(i) = fuel_tank_height - (fuelMassArr(i) / (fuel_density * ...
            densImp2Met) / tank_inner_area / 2);

        % Update propellant CoMs from nose
        loxCoM(i) = loxRelCoM(i) + heightToLox;
        fuelCoM(i) = fuelRelCoM(i) + heightToFuel;
        
        % Update total CoM from nose
        totCoM(i,2) = (loxCoM(i) * loxMassArr(i) + fuelCoM(i) * fuelMassArr(i) ...
                      + rocket_empty_CoM * mass_dry) / (loxMassArr(i) + fuelMassArr(i) ...
                      + mass_dry);
        totMass(i,2) = loxMassArr(i) + fuelMassArr(i) + mass_dry;

        % Update total simulation time
        totCoM(i,1) = tspan(i);
        totMass(i,1) = tspan(i);
    end
    % Update inertia
    com = totCoM(i,2);
    MoI(i,1,1) = 1/2 * totMass(i,2) * (tank_OD/2)^2;
    transverseMoI = 0;
    transverseMoI = transverseMoI + loxMassArr(i) * (1/4*(tank_OD/2)^2 + 1/12*lox_tank_height^2 + (loxCoM(i)-com)^2);
    transverseMoI = transverseMoI + fuelMassArr(i) * (1/4*(tank_OD/2)^2 + 1/12*fuel_tank_height^2 + (fuelCoM(i)-com)^2);
    transverseMoI = transverseMoI + engine_mass * (1/4*(tank_OD/2)^2 + 1/3*engine_height^2 + (engineHFore-com)^2);
    transverseMoI = transverseMoI + fin_can_mass * (1/4*(tank_OD/2)^2 + 1/3*fin_can_height^2 + (finCanHFore-com)^2);
    transverseMoI = transverseMoI + mid_mass * (1/4*(tank_OD/2)^2 + 1/3*mid_height^2 + (midAFHFore-com)^2);
    transverseMoI = transverseMoI + pressurant_tank_mass * (1/4*(tank_OD/2)^2 + 1/3*pressurant_tank_height^2 + (heHFore-com)^2);
    %the nose is currently a cylinder
    transverseMoI = transverseMoI + nose_mass * (1/4*(tank_OD/2)^2 + 1/3*nose_height^2 + (noseHFore-com)^2);
    
    MoI(i,2,2) = transverseMoI;
    MoI(i,3,3) = transverseMoI;
end

%% Optional CoM graphing output
if plot_CoM_graph == 1

    % Initialize arrays and text used for graphing
    totCoMTop = totCoM(:,2);
    loxRelCoMIndividual = lox_tank_height - loxRelCoM(:,1);
    fuelRelCoMIndividual = fuel_tank_height - fuelRelCoM(:,1);
    txtRock = ['Rocket height: ' num2str(total_height) ' meters'];
    txtLox = ['Lox tank height: ' num2str(lox_tank_height) ' meters'];
    txtFuel = ['Fuel tank height: ' num2str(fuel_tank_height) ' meters'];

    % Figure 1 - total CoM of the rocket
    figure(1)
    plot(tspan, totCoMTop, lineWidth=1.5)
    text(25,2.975,txtRock)
    grid on
    title("Center of Mass vs. Time")
    xlabel("Time [s]")
    ylabel("CoM from nose [m]")

    % Figure 2 - CoM of the lox tank
    figure(2)
    plot(tspan, loxRelCoM, lineWidth=1.5)
    hold on
    plot(tspan, loxRelCoMIndividual, lineWidth=1.5)
    text(25,0.45,txtLox)
    grid on
    title("loxCoM from top of tank")
    xlabel("Time [s]")
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
    xlabel("Time [s]")
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
