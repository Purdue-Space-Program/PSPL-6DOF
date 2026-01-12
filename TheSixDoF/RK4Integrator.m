function [out, mach, AoA, accel, cD, momentVector] = RK4Integrator(time, input, atmosphere, totCoM, totMass, InertMatrix, windData, rocket, settings, env)
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


pos = [input(1);input(2);input(3)]; % x y z ? - david

velocity = [input(4);input(5);input(6)];

omega = [input(7); input(8); input(9)];

quat = [input(10); input(11); input(12); input(13)];

components = values(rocket.ComponentList);

for idx = 1:length(components)
    if isa(components{idx}, 'PropulsionSystem')
        motor_fuck = components{idx};
    end
end


rocket_radius = rocket.OuterDiameter / 2;    % radius of rocket [m]
rocket_reference_area = pi*rocket_radius^2;                % reference area (m^2), as defined by RasAero (cross-sectional area)
thrust_magnitude = motor_fuck.Thrust;   % thrust of rocket in N.
body_vector = [1;0;0];       % vector in the body axis running through the nose.
chamber_exit_area = motor_fuck.ExitArea;     % exit area of the nozzle [m^2]
chamber_exit_pressure = motor_fuck.ExitPressure;      % exit pressure of the nozzle [Pa]
burn_time = motor_fuck.BurnTime;    
    

bodyVectorEarth = RotationMatrix(body_vector, quat, 1); % Body vector in inertial frame

%---------------- Get atmospheric Conditions (prefer Open-Meteo via env) ---
geoalt = real(pos(1));  % AGL [m]


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
windVector = windMag * [0;sin(windDir);cos(windDir)];


if settings.Wind == true
    freestream_velocity = velocity - windVector;
else
    freestream_velocity = velocity;
end

%---------------- Mass Update ----------------------------------------------

timeTableMass = totMass(:,1);
massTable = totMass(:,2);

[~,timeIndexMass] = min(abs(timeTableMass-time));
mass = massTable(timeIndexMass);
  
CoM_table_time_indexes = totCoM(:,1);
CoM_table = totCoM(:,2);

[~, timeIndexCoM] = min(abs(CoM_table_time_indexes-time));
CoM = CoM_table(timeIndexCoM);
% CoM = 1.5;

%---------------- Gravity force --------------------------------------------

if strcmpi(settings.Fidelity, "low")
    g = 9.8;
elseif strcmpi(settings.Fidelity,"medium") || strcmpi(settings.Fidelity,"high")
g = gravitywgs84(geoalt, env.LatLong(1), env.LatLong(2), 'Exact');
end

gravForce = mass * g * [-1;0;0];

% calculate the angle between the velocity vector and the rocket nose
AoA = acosd((dot(freestream_velocity,bodyVectorEarth)) / (norm(freestream_velocity) * norm(bodyVectorEarth)));
AoA(isnan(AoA)) = 0;
AoA = real(AoA);

% mach is used for airspeed dependent drag coefficient:
mach = norm(velocity) / a;

% read the coefficient of drag from RasAero data:
rasData = setAeroData(rocket);
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

pressure_thrust = thrust_magnitude + (chamber_exit_pressure - P) * chamber_exit_area;

if time <= burn_time
    thrustForceBody = pressure_thrust * body_vector;
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
    lift_coefficient = interp1([0,2],[cL0, cL2], AoA);
elseif AoA > 2 && AoA <= 4
    lift_coefficient = interp1([2,4],[cL2, cL4], AoA);
else
    slope = (cL4 - cL0) / 4;
    lift_coefficient = min(slope * AoA, 4);
end

% these act around the center of pressure, which is given in RasAero:
cPTable = rasData(1:300,7); % center of pressure in inches, defined from nose
cPTableMetric = cPTable / 39.37; %center of pressure in meters, defined from nose

[~, machIndex2] = min(abs(machTable-mach));

cP = cPTableMetric(machIndex2);

%find the magnitude of lift
lift = (0.5 * rho* norm(freestream_velocity)^2 * rocket_reference_area * lift_coefficient);

% do some vector math to find the lift direction:
liftDir = cross(cross(freestream_velocity, bodyVectorEarth), freestream_velocity) / norm(cross(cross(freestream_velocity,bodyVectorEarth),freestream_velocity));

liftForce = lift * liftDir;
liftForce(isnan(liftForce)) = 0;
liftForceBody = RotationMatrix(liftForce, quat, 0);


%---------------- Drag force -----------------------------------------------

% drag lies parellel and opposite to the velocity vector

%determine the direction and magnitude of the drag force
drag_direction = -freestream_velocity / norm(freestream_velocity);
% implement a simple drag polar model for drag increase with AoA:
cD = cD + 0.1*(lift_coefficient)^2;

drag_magnitude = (0.5 * rho * norm(freestream_velocity)^2 * rocket_reference_area * cD);
drag_force = drag_direction * drag_magnitude;
drag_force(isnan(drag_force)) = 0; % turn all NaN values to zero
dragForceBody = RotationMatrix(drag_force, quat, 0);


%---------------- Parachute ------------------------------------------------

persistent tDrogueOffset

if velocity(1) < 0

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
        calculateParachuteForce(pos, freestream_velocity, rho, env, tDrogue);

else
    paraDragForce = [0;0;0];
end

% Output magnitude of parachute force for reefing study
%fprintf("Para drag force is: %.2f\n", norm(paraDragForce));

% Convert to body frame for moments calculations
paraDragForceBody = RotationMatrix(paraDragForce, quat, 0);


%---------------- Total Forces ---------------------------------------------

forceVector = gravForce + thrustForceEarth + drag_force + liftForce + paraDragForce;

accel = forceVector / mass;

%---------------- Stability Caliber ----------------------------------------

% difference between CoM and cP divided by diameter of the rocket
fprintf("\nCoM: %.3f, cP: %.3f", CoM, cP)
fprintf("\nStability caliber: %.3f\n", (cP - CoM) / rocket.OuterDiameter);

%---------------- Moments --------------------------------------------------

% pull the moments from the CoM MoI data:

Ixx = InertMatrix(timeIndexMass,1,1);
Iyy = InertMatrix(timeIndexMass,2,2);
Izz = InertMatrix(timeIndexMass,3,3);

I = [Ixx, 0, 0;
     0, Iyy, 0;
     0, 0, Izz];

AeroMomentArm = (CoM - cP) * body_vector; % define the length of the moment arm in the body frame
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
coefficientLift = 5e-6 * missAlpha * 0;

forceRoll = 3 / 2 * coefficientLift * rho * norm(velocity)^2;
rollMomentBody = (rocket_radius + finCpLocation) * forceRoll * body_vector;

momentVector = liftMomentBody + dragMomentBody + rollMomentBody + paraMomentBody;

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

out = [velocity;accel;alpha;quatRates];
end