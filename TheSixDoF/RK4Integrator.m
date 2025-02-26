function [out, mach, AoA, accel] = RK4Integrator(time, input, aeroArray, rocketParams, atmosphere, totCoM, totMass, J, wind, params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% aeroData - Array of mach, Cd, and Cp from CdModel and CpModel
% totCoM - Array of center of mass measurements at different time steps
%          [s|m]
% totMass - Array of total mass values at different time steps [s|kg]
% J - Moment of Inertia of the rocket [m^4]
% wind - Array of data with wind information
% params - extraneous parameters to be passed into function
%
% Outputs:
% out = derivative of state vector [m/s|m/s^2|rad/s^2|unitless^2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos = [input(1);input(2);input(3)];

vel = [input(4);input(5);input(6)];

omega = [input(7); input(8); input(9)];

quat = [input(10); input(11); input(12); input(13)];

radius = rocketParams(2,3) / 2;    % radius of rocket [m]
A = (pi * radius ^ 2) + (rocketParams(2,5) * rocketParams(2,6) * rocketParams(2,11));          % reference area (m^2)
thrustMag = rocketParams(2,14);  % average thrust of rocket [N].
burnTime = rocketParams(2,15);        % burn time of motor [s]
bodyVector = [1;0;0]; % vector in the body axis running through the nose.
ExitA = rocketParams(2,16);    % exit area of the nozzle [m^2]
ExitP = rocketParams(2,17);      % exit pressure of the nozzle [Pa]

if nargin == 6
    thrustMag = params(1);
end


bodyVectorEarth = RotationMatrix(bodyVector, quat, 1); % Body vector in inertial frame

%% atmospheric conditions:
height = pos(1);

% long monte carlo can result in complex numbers here, take real component
height = real(height);

% get atmospheric parameters
%[~, a, P, rho] = atmosisa(height);
heightIndex = int32(height)+1;
a = atmosphere(heightIndex, 1);
rho = atmosphere(heightIndex,2);
P = atmosphere(heightIndex, 3);

%% Wind:
windAlt = wind(:,1);
windMagList = wind(:,2);
windDirList = wind(:,3);

%% Mass Update:
timeTableMass = totMass(:,1);
massTable = totMass(:,2);

[~,timeIndexMass] = min(abs(timeTableMass-time));
mass = massTable(timeIndexMass);

%% Wind calcs:
[~, heightIndex] = min(abs(windAlt-height));

windDir = windDirList(heightIndex);
windMag = windMagList(heightIndex);
windVector = windMag * [sin(windDir);cos(windDir);0];

windVel = vel - windVector;
  
%% Center of mass update
timeTableCoM = totCoM(:,1);
CoMTable = totCoM(:,2);

[~, timeIndexCoM] = min(abs(timeTableCoM-time));
CoM = CoMTable(timeIndexCoM);

%% Gravitational Force:
gravForce = mass * constant.g * [-1;0;0];

%% Thrust Forces:
% thrust lies along long axis of the rocket [1;0;0], which we then convert into
% earth frame

presThrust = thrustMag + (ExitP - P) * ExitA;

if time <= burnTime
    thrustForceBody = presThrust * bodyVector;
    thrustForceEarth = RotationMatrix(thrustForceBody, quat, 1);
else
    thrustForceEarth = [0;0;0];
end

%% Drag Forces:
% drag lies parellel and opposite to the velocity vector

% calculate the angle between the velocity vector and the rocket nose
AoA = acosd((dot(windVel,bodyVectorEarth)) / (norm(windVel) * norm(bodyVectorEarth)));
AoA(isnan(AoA)) = 0;

% mach is used for airspeed dependent drag coefficient:
mach = norm(vel) / a;

% read the coefficient of drag from RasAero data:
machTable = aeroData(1:300,1); % mach values from 0.01 to 3
cDTable = aeroData(1:300,3); % coefficient of drag

% find cD matching the closest mach value to table
[~, machIndex] = min(abs(machTable-mach));
cD = cDTable(machIndex);

%determine the direction and magnitude of the drag force
dragDir = -windVel / norm(windVel);
dragMag = (0.5 * rho * norm(windVel)^2 * A * cD);

% implement a simple quadratic model for drag increase with AoA:
dragMag = min(dragMag + 0.075*dragMag*(AoA)^2, 3 * dragMag);

dragForce = dragDir * dragMag;
dragForce(isnan(dragForce)) = 0;
dragForceBody = RotationMatrix(dragForce, quat, 0);


%% Parachute Drag

% Other Constants
deformVal = 0.70;                   % Deformation value of inflated chute area
                                    % *This is an approximate value due to
                                    % the difficulty to calculate it*

% Parachute parameters
drogueCd = 0.97;                    % cD for the drogue chute
drogueDia = (25/6) * constant.ft2m; % drogue  diameter [m]

mainCd = 2.2;                       % cD for the main parachute
mainDia = (97/6) * constant.ft2m;   % main chute diameter [m]
mainDeployAlt = 304.8;              % main chute deployment altitude [m]

if vel(1) < 0
    % Drogue deployment only
    if pos(1) > mainDeployAlt
        totDia = drogueDia;
        totCd = drogueCd;

    % Main deployment
    else        
        % Assume instantaneous opening for testing purposes
        totDia = (drogueDia + mainDia);
        totCd = drogueCd + mainCd;
    end   

    vertArea = deformVal * pi * (0.5 * totDia) ^ 2;
    forceMag = 0.5 * totCd * vertArea * rho * norm(vel) ^ 2;
    paraDragForce = forceMag .* (-vel ./ norm(vel)); 

else
    paraDragForce = [0;0;0];
end

% Output magnitude of parachute force for reefing study
%fprintf("Para drag force is: %.2f\n", norm(paraDragForce));

% Convert to body frame for moments calculations
paraDragForceBody = RotationMatrix(paraDragForce, quat, 0);


%% Lift Forces:
% lift forces lie perpendicular to the velocity, these are the most
% difficult to calculate accurately

% a simple linear model for the coefficient of lift wrt on AoA is
% being used right now, in the future VSPaero or CFD data may be used.
% These values are loosely based on DATCOM / RasAero data we have
% previously gathered.

cL = min(1/8 * AoA, 2);

% these act around the center of pressure, which is given in RasAero,
cPTable = aeroData(1:300,13); % center of pressure in inches, defined from nose
cPTableMetric = cPTable / 39.37; %center of pressure in meters, defined from nose

[~, machIndex2] = min(abs(machTable-mach));

cP = cPTableMetric(machIndex2);

% do some vector math to find the lift direction:
lift = (0.5 * rho* norm(windVel)^2 * A * cL);

liftDir = cross(cross(windVel, bodyVectorEarth), windVel) / norm(cross(cross(windVel,bodyVectorEarth),windVel));

liftForce = lift * liftDir;
liftForce(isnan(liftForce)) = 0;
liftForceBody = RotationMatrix(liftForce, quat, 0);

%% Total Forces:
forceVector = gravForce + thrustForceEarth + dragForce + liftForce + paraDragForce;
%accel:
accel = forceVector / mass;

%% Stability Caliber Calculations
% difference between CoM and cP divided by diameter of the rocket
%fprintf("Stability caliber: %.3f\n", abs(CoM - cP) / 0.168275);

%% Moments:
% just putting in values based on a cylinder moment of inertia for now.
Jxx = 1;
Jyy = 200;
Jzz = 200;

%% Aerodynamic Moments:

AeroMomentArm = (CoM - cP) * bodyVector; % define the length of the moment arm in the body frame
%ParaMomentArm = CoM * bodyVector'; % define the length of the moment arm of the parachute in the body frame

% lift:
%liftForce = RotationMatrix(liftForce, angle, 0); % define the lift force in body axes
liftMomentBody = cross(AeroMomentArm,liftForceBody);
% drag:
dragMomentBody = cross(AeroMomentArm,dragForceBody);

% Needs work to be correct
% paraMomentBody = cross(ParaMomentArm, paraDragForceBody);
% AoAPara = acosd((dot(-vel,bodyVectorEarth)) / (norm(vel) * norm(bodyVectorEarth)));
% AoAPara(isnan(AoAPara)) = 0;
% 
% paraMomentBody = exp(0.1*(38-time)) * 0.15 .* paraMomentBody + 0.3 * -sind(AoAPara) .* paraMomentBody;


%% Roll Moment Test:

finCpLocation = 0.02486256; % 1/3 of the span of fins [m]
missAlpha = 1; % [degrees]
coefficiantLift = 0.0005 * missAlpha;

missAlpha = 0; %degrees
coefficiantLift = 0.0005*missAlpha;

forceRoll = 3 / 2 * coefficiantLift * rho * norm(vel)^2;
rollMomentBody = (radius + finCpLocation) * forceRoll * bodyVector;

momentVector = liftMomentBody + dragMomentBody + rollMomentBody; %+ paraMomentBody;

% use euler equations to find the final moments:

momentVector(1) = momentVector(1) - omega(2)*omega(3)*(Jzz-Jyy);
momentVector(2) = momentVector(2) - omega(1)*omega(3)*(Jxx-Jzz);
momentVector(3) = momentVector(3) - omega(1)*omega(2)*(Jyy-Jxx);

alpha = inv(J) * momentVector;
alpha(isnan(alpha)) = 0;

wx = omega(1);
wy = omega(2);
wz = omega(3);

B = [0, -wx, -wy, -wz;
     wx, 0, wz, -wy;
     wy, -wz, 0, wx;
     wz, wy, -wx, 0];

quatRates = 0.5 * B * quat;

out = [vel;accel;alpha;quatRates];



end