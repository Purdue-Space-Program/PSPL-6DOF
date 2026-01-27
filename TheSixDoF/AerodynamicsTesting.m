% load the rasAero data

rasAeroData = readmatrix("Inputs/RASAero/RasAeroDataCulled2.csv");

% get the parameters for zero angle of attack

rasMach = rasAeroData(1:300,1);

rasCd = rasAeroData(1:300,3);

rasCl = rasAeroData(1:300,6);

rasCp = rasAeroData(1:300,7); %this is in inches, from the nose

rasCp = rasCp/12 * constant.FT_TO_M; % convert from inches to meters

% load the rocket object and compute its aerodynamics:
rocket = load("Inputs\Saved Rockets\CMS.mat");

rocket = rocket.rocketObj;

% input the rocket into the functions and calculate aerodynamic
% characteristics:

drag = computeDrag(rocket, mach, alpha);
lift = computeLift(rocket, mach, alpha);
% cp = computeCp(rocket, mach, alpha);


% plot the drag as a function of mach number:
figure;
plot(rasMach, rasCd);
xlabel('Mach Number');
ylabel('Drag Coefficient ($C_D$)');
title('Drag Coefficient vs Mach Number');
grid on;

% plot the lift as a function of mach number:
figure;
plot(rasMach, rasCl);
xlabel('Mach Number');
ylabel('Drag Coefficient ($C_L$)');
title('Drag Coefficient vs Mach Number');
grid on;

% plot the center of pressure location as a function of mach number:
figure;
plot(rasMach, rasCp);
xlabel('Mach Number');
ylabel('Center of Pressure Location (m)');
title('Center of Pressure vs Mach Number');
grid on;
