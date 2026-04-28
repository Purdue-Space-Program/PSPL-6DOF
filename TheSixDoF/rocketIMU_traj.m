%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSP FLIGHT DYNAMICS:
%
% Title: Monte
% Author: Hudson Reynolds
%
% Description: Runs the 6DoF with simulated inertial navigation
%
% Inputs: N/A
%a
% Outputs:
% see subfunctions for specific outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
close all

% add all the folders
addpath(genpath('Classes'))

% % Make a GPS with measurement update:
gps = GNSS("GPS",2, [4,4,25], .1, 0);

% -- set up the simulation ---:
env = Environment();
settings = IntegratorSettings("apogee", 0.01, "high");
env = getLocalWeather(env);

%% Run the sim:

filename = ['Inputs',filesep,'Saved Rockets',filesep,'CMS.mat'];

rocket = load(filename);
rocket = rocket.rocketObj;

% create a gyroscope and accelerometer:
accel = rocket.ComponentList.values{6};
gyro = rocket.ComponentList.values{7};

gyro.SamplingRate = 1/100;
gyro.Variance = .0035*pi/180*ones(3,1); % from VN-200
gyro.Bias0 = 0;
gyro.ScaleFactor = 1.001;

accel.SamplingRate = 1/100;
accel.Variance = 1e-3*ones(3,1);
accel.Bias0 = 1e-4;
accel.ScaleFactor = 1.001;

% add the gyro and accel to the rocket
rocket.modifyComponent(gyro)
rocket.modifyComponent(accel)

components = values(rocket.ComponentList);


time = 200;

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

tic;
[timeArray, out] = ode45(@(time,input) RK4Integrator(time,input,atmosphere, ...
totCoM,totMass,MoI,windData, rocket, settings, env), tspan, Init, opt);
toc;


%% Outputs:

% create a struct which contains all of the output information:
outputStruct = struct;
outputStruct.time = timeArray;

% output additional arrays from the integrator
for k = 1:numel(timeArray)
    [~, outputStruct.mach(k,1), outputStruct.AoA(k,1), outputStruct.accel(k,:), ...
        outputStruct.accelBody(k,:), outputStruct.cD(k,:), moment(k,:), outputStruct.g(k,:)] = ...
        RK4Integrator(timeArray(k), out(k,:), ...
        atmosphere,totCoM, totMass, MoI, windData, rocket, settings, env);
end


if settings.Outputs == true

    % parse rk4 outputs:
    posArray = out(:,1:3);
    velArray = out(:,4:6);
    omega = out(:,7:9);
    quatArray = out(:,10:13);

    outputStruct.pos = posArray;
    outputStruct.vel = velArray;
    outputStruct.omega = omega;
    outputStruct.quat = quatArray;

    E = wgs84Ellipsoid;
    [lats,longs, ~] = enu2geodetic(out(:,1),out(:,2),out(:,3),env.LatLong(1),env.LatLong(2),E.SemimajorAxis,E);


    uif = uifigure;
    g = geoglobe(uif);
    
    geoplot3(g, lats, longs, out(:,3), 'r-', LineWidth= 1)
    campos(g,env.LatLong(1)-0.1,env.LatLong(2)-0.1,15000)
    campitch(g,-30)
    camheading(g,45)
    

    % find end conditions for graphs / animations
    endTime = length(outputStruct.AoA) * settings.Timestep;

    % grab parameters at max Q and off the rail
    [maxVel, maxIndex] = max(out(:,4));
    maxqAccel = outputStruct.accel(maxIndex,1);
    maxqpos = posArray(maxIndex,1);

    machTable = rasData(1:300,1);
    cdTable = rasData(1:300,3);
    maxqMach = outputStruct.mach(maxIndex);
    [~, maxqMachIndex] = min(abs(machTable-maxqMach));
    maxqCD = cdTable(maxqMachIndex);

    % [~, railIndex] = min(abs(posArray(1:100,1)-env.railHeight));
    % railVel = out(railIndex,4);
    % railAccel = outputStruct.accel(railIndex,1);
    % 
    % railMach = outputStruct.mach(railIndex);
    % [~, railMachIndex] = min(abs(machTable-railMach));
    % railCD = cdTable(railMachIndex);

    apogee = max(posArray(:,3));

    % fprintf("\nParameters at Max Q:\n")
    % fprintf(" Velocity: %.2f m/s\n Mach: %.3f\n Acceleration: %.3f m/s^2\n Drag Coefficient: %.4f\n",maxVel, maxqMach, maxqAccel, maxqCD);
    % fprintf("Off-Rail Parameters:\n")
    % fprintf(" Velocity: %.2f m/s\n Mach: %.3f\n Acceleration: %.3f m/s^2\n Drag Coefficient: %.4f\n",railVel, railMach, railAccel, railCD);
    % fprintf("Rocket Apogee (AMSL): %.2f m\n", apogee)
    % fprintf("Rocket Apogee (AGL): %.2f m\n", apogee - env.Elevation)

    

    % Plotting:
    gyroObj = [];
    for i = 1:numel(components)
        c = components{i};
        if isa(c, 'Gyroscope')
            gyroObj = c;
            break
        end
    end


    % Gyroscope plot
    if ~isempty(gyroObj)
        omegaTrue = outputStruct.omega;
        [omegaIMU,omegaIMUTime] = gyroObj.GyroscopeMeasurement(omegaTrue, timeArray(end),settings.Timestep);
        %gyroObj.plotMeasurementHistory(timeArray, omegaTrue, omegaMeas);
    end
    
    accelObj = [];
    for i = 1:numel(components)
        c = components{i};
        if isa(c, 'Accelerometer')
            accelObj = c;
            break
        end
    end

    % Accelerometer plot
    if ~isempty(accelObj)
        accelTrue = outputStruct.accelBody;
        [accelIMU, accelIMUTime] = accelObj.AccelerometerMeasurement(accelTrue, timeArray(end),settings.Timestep);
        %accelObj.plotMeasurementHistory(timeArray, accelTrue, accelMeas)
    end

    % get the GPS pings:
    posTrue = outputStruct.pos;
    [GPSpos, GPSposTime] = gps.GNSSMeasurement(posTrue, timeArray(end),settings.Timestep);

    % Rocket Trajectory Plot:
    figure;
    plot3(posArray(1:int32(endTime / settings.Timestep),1), posArray(1:int32(endTime / settings.Timestep),2), posArray(1:int32(endTime / settings.Timestep),3))
    hold on
    plot3(GPSpos(:,1), GPSpos(:,2), GPSpos(:,3), 'go')
    % plot3(posArray(1:endTime / dt,3), posArray(1:endTime / dt,2), zeros(endTime / dt), '--')
    % plot3(posArray(1:endTime / dt,3), zeros(endTime / dt), posArray(1:endTime / dt,1), '--')
    % plot3(zeros(endTime / dt), posArray(1:endTime / dt,2), posArray(1:endTime / dt,1), '--')
    view(43,24);
    xlabel('Dist North (m)');
    ylabel('Dist East (m)');
    zlabel('Height (m)');
    %axis equal;
    grid minor;

end

IMU_data = struct;

IMU_data.accelometer = accelIMU;
IMU_data.accel_time = accelIMUTime;
IMU_data.gyroscope = omegaIMU;
IMU_data.gyro_time = omegaIMUTime;
IMU_data.GPS = GPSpos;
IMU_data.GPS_time = GPSposTime;

IMU_data.init_cond = Init([1:6,10:13]); %pos, vel, quat
IMU_data.ref_traj = posArray;
IMU_data.ref_time = timeArray;
IMU_data.grav = outputStruct.g;

IMU_data.accel_info = accel;
IMU_data.gyro_info = gyro;
IMU_data.gps_info = gps;

save("Rocket_IMU_data.mat", "IMU_data")
