function [out, mach, AoA, accel, specific_force_body, cD, momentVector, g] = RK4Integrator(time, input, atmosphere, totCoM, totMass, InertMatrix, windData, rocket, settings, env)
% PSP FLIGHT DYNAMICS:
%
% Title: RK4Integrator
% Author: Hudson Reynolds - Created: 9/21/2024
% Editors: Preston Wright (implemented variable CoM/stability caliber)
%          - 9/28/2024
%
% Description: This is the integration function to be used in ode45. This
% computes all funciton derivatives and differential equations for the
% translational and rotational dynamics.
%
% Inputs: 
% time - current simulation time [s]
% input - Array of position, velocity, rotational velocity, and quaternions
%         [m|m/s|rad/s|unitless] 
% rasData - Array of data pulled from ras aero [lots of units]
% totCoM - Array of center of mass measurements at different time steps
%          [s|m]
% totMass - Array of total mass values at different time steps [s|kg]
% J - Moment of Inertia of the rocket [m^4]
% wind - Array of data with wind information
% params - extraneous parameters to be passed into function
%
% Outputs:
% out = derivative of state vector [m/s|m/s^2|rad/s^2|unitless^2]

%---------------- Parse Function Inputs -------------------------------------


pos = [input(1);input(2);input(3)];

vel = [input(4);input(5);input(6)];

omega = [input(7); input(8); input(9)];

quat = [input(10); input(11); input(12); input(13)];

components = values(rocket.ComponentList);

for idx = 1:length(components)
    if isa(components{idx}, 'PropulsionSystem')
        motor = components{idx};
    end
end

radius = rocket.OuterDiameter / 2;    % radius of rocket [m]
A = pi*radius^2;                % reference area (m^2), as defined by RasAero (cross-sectional area)
thrustMag = motor.Thrust;   % thrust of rocket in N.
bodyVector = [1;0;0];       % vector in the body axis running through the nose.
ExitA = motor.ExitArea;     % exit area of the nozzle [m^2]
ExitP = motor.ExitPressure;      % exit pressure of the nozzle [Pa]
burnTime = motor.BurnTime;


bodyVectorEarth = RotationMatrix(bodyVector, quat, 1); % Body vector in inertial frame

%---------------- Get atmospheric Conditions (prefer Open-Meteo via env) ---
geoalt = real(pos(3));  % AGL [m]


% get the atmosphere data (update to get the environment):

% atmosphere (height, temp, a, pres, rho)

[~, atmosIdx] = min(abs(env.Atmosphere(:,1) - geoalt));

rho = env.Atmosphere(atmosIdx,5);
a = env.Atmosphere(atmosIdx,3);
T = env.Atmosphere(atmosIdx,2);
P = env.Atmosphere(atmosIdx,4);

%---------------- Wind -----------------------------------------------------

windAlt = windData(:,1);
windMagList = windData(:,2);
windDirList = windData(:,3);

[~, heightIndex] = min(abs(windAlt-geoalt));

windDir = windDirList(heightIndex);
windMag = windMagList(heightIndex);
windVector = windMag * [sin(windDir);cos(windDir);0];


if settings.Wind == true
    freestreamVel = vel - windVector;
else
    freestreamVel = vel;
end

%---------------- Mass Update ----------------------------------------------

timeTableMass = totMass(:,1);
massTable = totMass(:,2);

[~,timeIndexMass] = min(abs(timeTableMass-time));
mass = massTable(timeIndexMass);
  
timeTableCoM = totCoM(:,1);
CoMTable = totCoM(:,2);

[~, timeIndexCoM] = min(abs(timeTableCoM-time));
CoM = CoMTable(timeIndexCoM);

%---------------- Gravity force --------------------------------------------

if strcmpi(settings.Fidelity, "low")
    g = 9.8;
elseif strcmpi(settings.Fidelity,"medium") || strcmpi(settings.Fidelity,"high")
g = gravitywgs84(geoalt, env.LatLong(1), env.LatLong(2), 'Exact');
end

gravForce = mass * g * [0;0;-1];

% calculate the angle between the velocity vector and the rocket nose
AoA = acosd((dot(freestreamVel,bodyVectorEarth)) / (norm(freestreamVel) * norm(bodyVectorEarth)));
AoA(isnan(AoA)) = 0;
AoA = real(AoA);

% mach is used for airspeed dependent drag coefficient:
mach = norm(vel) / a;

% read the coefficient of drag from RasAero data:
rasData = rocket.AeroData;
machTable = rasData(1:300,1); % mach values from 0.01 to 3
cDTable = rasData(1:300,3); % coefficient of drag

% cL vs. AoA:
cLAoA0 = rasData(1:300, 6);
cLAoA2 = rasData(301:600, 6);
cLAoA4 = rasData(601:900, 6);

% find cD matching the closest mach value to table
[~, machIndex] = min(abs(machTable-mach));
cD = cDTable(machIndex);

%---------------- Thrust force ---------------------------------------------

% thrust lies along long axis of the rocket [1;0;0], which we then convert into
% earth frame

presThrust = thrustMag + (ExitP - P) * ExitA;

if time <= burnTime
    thrustForceBody = presThrust * bodyVector;
    thrustForceEarth = RotationMatrix(thrustForceBody, quat, 1);
else
    thrustForceEarth = [0;0;0];
end

%---------------- Lift force -----------------------------------------------

% lift forces lie perpendicular to the velocity, these are the most
% difficult to calculate accurately

% a simple linear model for the coefficient of lift wrt on AoA is
% being used right now, in the future VSPaero or CFD data may be used.
% These values are loosely based on DATCOM / RasAero data we have
% previously gathered. Currently there is no dependence on mach either,
% which should be implemented at some point

% find the lift coeff. at closest mach number
cL0 = cLAoA0(machIndex);
cL2 = cLAoA2(machIndex);
cL4 = cLAoA4(machIndex);

if AoA <= 2
    cL = interp1([0,2],[cL0, cL2], AoA);
elseif AoA > 2 && AoA <= 4
    cL = interp1([2,4],[cL2, cL4], AoA);
else
    slope = (cL4 - cL0) / 4;
    cL = min(slope * AoA, 4);
end

% these act around the center of pressure, which is given in RasAero:
cPTable = rasData(1:300,7); % center of pressure in inches, defined from nose
cPTableMetric = cPTable / 39.37; %center of pressure in meters, defined from nose

[~, machIndex2] = min(abs(machTable-mach));

cP = cPTableMetric(machIndex2);

%find the magnitude of lift
lift = (0.5 * rho* norm(freestreamVel)^2 * A * cL);

% do some vector math to find the lift direction:
liftDir = cross(cross(freestreamVel, bodyVectorEarth), freestreamVel) / norm(cross(cross(freestreamVel,bodyVectorEarth),freestreamVel));

liftForce = lift * liftDir;
liftForce(isnan(liftForce)) = 0;
liftForceBody = RotationMatrix(liftForce, quat, 0);


%---------------- Drag force -----------------------------------------------

% drag lies parellel and opposite to the velocity vector

%determine the direction and magnitude of the drag force
dragDir = -freestreamVel / norm(freestreamVel);
% implement a simple drag polar model for drag increase with AoA:
cD = cD + 0.1*(cL)^2;

dragMag = (0.5 * rho * norm(freestreamVel)^2 * A * cD);
dragForce = dragDir * dragMag;
dragForce(isnan(dragForce)) = 0;
dragForceBody = RotationMatrix(dragForce, quat, 0);


%---------------- Parachute ------------------------------------------------

persistent tDrogueOffset

if vel(3) < 0

    % initialize a new time variable for the drogue deployment which
    % increases the length of the drogue as a function of time.
    if isempty(tDrogueOffset)
        tDrogueOffset = time;
        tDrogue = 0;
    else
        tDrogue = time - tDrogueOffset;
        if tDrogue < 0
            tDrogue = 0;
        end
    end
    
    paraDragForce = ...
        calculateParachuteForce(pos, freestreamVel, rho, env, tDrogue);

else
    paraDragForce = [0;0;0];
end

% Output magnitude of parachute force for reefing study
%fprintf("Para drag force is: %.2f\n", norm(paraDragForce));

% Convert to body frame for moments calculations
paraDragForceBody = RotationMatrix(paraDragForce, quat, 0);


%---------------- Total Forces ---------------------------------------------

forceVector = gravForce + thrustForceEarth + dragForce + liftForce + paraDragForce;
%accel:
accel = forceVector / mass;

g = gravForce/mass;

specific_force_world = accel - g;

specific_force_body = RotationMatrix(specific_force_world, quat, 0);

%---------------- Stability Caliber ----------------------------------------

% difference between CoM and cP divided by diameter of the rocket
%fprintf("Stability caliber: %.3f\n", abs(CoM - cP) / rocket.OuterDiameter);

%---------------- Moments --------------------------------------------------

% pull the moments from the CoM MoI data:

Ixx = InertMatrix(timeIndexMass,2);
Iyy = InertMatrix(timeIndexMass,3);
Izz = InertMatrix(timeIndexMass,4);

I = [Ixx, 0, 0;
     0, Iyy, 0;
     0, 0, Izz];

AeroMomentArm = (CoM - cP) * bodyVector; % define the length of the moment arm in the body frame
%ParaMomentArm = CoM * bodyVector'; % define the length of the moment arm of the parachute in the body frame

liftMomentBody = cross(AeroMomentArm,liftForceBody);
dragMomentBody = cross(AeroMomentArm,dragForceBody);

% This parachute moment should come out of the nose of the rocket. Update
% this to include this behavior and introduce stabilization in post-apogee.

paraMomentArm = [CoM;0;0];

paraMomentBody = cross(paraMomentArm,paraDragForceBody);

%---------------- Roll Moments ---------------------------------------------

finCpLocation = 0.02486256; % 1/3 of the span of fins [m]
missAlpha = 0.1; % [degrees]
coefficientLift = 5e-6 * missAlpha;

forceRoll = 3 / 2 * coefficientLift * rho * norm(vel)^2;
rollMomentBody = (radius + finCpLocation) * forceRoll * bodyVector;

momentVector = liftMomentBody + dragMomentBody + rollMomentBody; %+ paraMomentBody;

% use euler equations to find the final moments:

omegaX = omega(1);
omegaY = omega(2);
omegaZ = omega(3);

momentVector(1) = momentVector(1) - omegaY*omegaZ*(Izz-Iyy);
momentVector(2) = momentVector(2) - omegaX*omegaZ*(Ixx-Izz);
momentVector(3) = momentVector(3) - omegaX*omegaY*(Iyy-Ixx);

alpha = inv(I) * momentVector;
alpha(isnan(alpha)) = 0;

omegaX = omega(1);
omegaY = omega(2);
omegaZ = omega(3);

B = [0, -omegaX, -omegaY, -omegaZ;
     omegaX, 0, omegaZ, -omegaY;
     omegaY, -omegaZ, 0, omegaX;
     omegaZ, omegaY, -omegaX, 0];

quatRates = 0.5 * B * quat;

out = [vel;accel;alpha;quatRates];
end