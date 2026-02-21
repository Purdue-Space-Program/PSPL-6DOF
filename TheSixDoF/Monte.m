%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSP FLIGHT DYNAMICS:
%
% Title: Monte
% Author: Hudson Reynolds
%
% Description: Runs the 6DoF with monte carlo for landing location
%
% Inputs: N/A
%
% Outputs:
% see subfunctions for specific outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all
TRIALS = 50;

% settings:

% run Mojave:
env = Environment();

settings = IntegratorSettings("full", 0.1, "low");

env = getLocalWeather(env);

%% Run the sim:

rocket = load('Inputs\Saved Rockets\CMS.mat');

rocket = rocket.rocketObj;

components = values(rocket.ComponentList);

if strcmpi('burnout', settings.EndCondition)
    time = rocket.BurnTime;
elseif ~isnan(str2double(settings.EndCondition))
    time = round(str2double(settings.EndCondition),1);
else
    time = 400;
end

arrayLength = (time / settings.Timestep);
tspan = linspace(0,time,arrayLength+1);

% set the initial position in ENU frame(x,y,z). Accounts for starting elevation.
pos = [0;0;env.Elevation];

% set the initial velocity in ENU frame (xdot,ydot,zdot).
vel = [0;0;0];

% initial angle (z angle, y angle, x angle) - following 3-2-1 sequence
angleVector = [0;-pi/2;0];

% initial rotation rate (x rate, y rate, z rate)
omega = [0;0;0];

% initalize the quaternion based on the euler angle input:
quatVector = eul2quat(angleVector.', "ZYX").';

% initial state vector
Init = [pos;vel;omega;quatVector];

% import aerodynamics data
rasData = rocket.AeroData;

% import wind data (prefer Open-Meteo via env, fallback to parser)
if (isstruct(env) && isfield(env,'WindData')) || (isobject(env) && isprop(env,'Wind'))
    windData = env.Wind;   % [alt_m, speed_mps, dir_rad]
else
    windData = wind.parseWind();
end


% import atmosphere;
atmosphere = env.Atmosphere;

% create an array of the center of mass, mass, and moment of inertia of the
% rocket
[totMass, totCoM, MoI] = VariableCoMMoI(rocket);

% additional options for RK4 (stop after reaching final condition)
opt = odeset('Events', @(tspan, Init) stoppingCondition(tspan, Init, settings.EndCondition, env), ...
    'RelTol', settings.relTol, 'AbsTol', settings.absTol);

% nominal thrust

thrustNom = 3.7929e+03;
burnNom = 13;


% do monte carlo, vary the launch angle and wind for number of trials
for trial = 1:TRIALS
    % randomize the angle vector by 
    zAng = 2*pi*rand(1);
    yAng = -pi/2 + 2*0.0175*randn(1);
    xAng = 0;


    launchAngle = [zAng; yAng; xAng];
    quatVector = eul2quat(launchAngle.', "ZYX").';
    Init = [pos;vel;omega;quatVector];
    % Run the simulation for the current trial

    % vary the thrust of the rocket
    components = rocket.ComponentList.values;

    prop = components{5};

    thrustVar = 100*randn(1);

    prop.Thrust = thrustNom + thrustVar;
    prop.BurnTime = prop.BurnTime - thrustVar/291.76;

    % update the rocketObj
    rocket.modifyComponent(prop)
    
    fieldName = sprintf('test%d', trial);


    [timeArray, out.(fieldName)] = ode45(@(time,input) RK4Integrator(time,input,atmosphere, ...
    totCoM,totMass,MoI,windData, rocket, settings, env), tspan, Init, opt);

    disp(trial)
end




%% Outputs:

wgs84 = wgs84Ellipsoid;

for idx = 1:TRIALS

    fieldName = sprintf('test%d', idx);

    result = out.(fieldName);

    % get the max apogee for that trial:
     apogee(idx) = max(result(:,3));

     % get the ending location east and north:
     northLocation(idx) = result(end, 2);
     eastLocation(idx) = result(end, 1);

     % convert to location and geoplot the final landing spot

     % uif = uifigure;
     % g = geoplot(uif);


     % use the starting coords and displacement for final coords:
     [lat(idx), long(idx), h] = enu2geodetic(eastLocation(idx),northLocation(idx), env.Elevation, env.LatLong(1), env.LatLong(2), env.Elevation, wgs84);

end

     % plot the final lat long on a geoplot

     g = geoaxes();
     set(g, 'Basemap', 'satellite')
     geoplot(g, lat, long, 'rx')
     hold on
     geoplot(g, env.LatLong(1), env.LatLong(2), 'bo')
     title('Monte Carlo with $T=\mathcal{N}(3792.9,100)$ and 1 deg launch angle randomization')

     % output the apogees:
    disp(mean(apogee))


