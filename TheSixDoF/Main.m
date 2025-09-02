% PSP FLIGHT DYNAMICS:
%
% Title: MainRK4
% Author: Hudson Reynolds - Created: 9/21/2024
% Last Modified: 6-20-2025
%
% Description: This is the overarching function that runs the 6-DoF,
% calling all neccesary functions to run the simulation. The overarching
% simulation structure uses an RK4 structure using ODE45.
%
% Inputs: Changes to the simulation are made by changing simulation
% settings
%
% Outputs:
% Graph, value, and file outputs. See subfunctions for specific outputs.


% ---------------- KNOWN ISSUES (WIP) ----------------------------------------

% post apogee attitude dynamics are not fully finished, attitude in this
% regime is likely incorrect
%
% non-zero AoA calculations for lift and drag need better modelling
%
% Monte Script is not up to date with main.


%---------------- Initialization ---------------------------------------------

% clear the console and figures before running the code:
clear;clc;close all force

% Create a rocket object, the default values are for CMS:
rocket = Rocket.Rocket();

%---------------- Sensor Definition ------------------------------------------

% Make a really bad altimeter for testing
altimeter = Sensor.Altimeter("Altimeter", 0.25, 20^2,.5, 5, .01);

% Make a GPS with measurement update:
gps = Sensor.GNSS("GPS",2, 3^2, .1, 0);

% Make a magnetometer:
mag = Sensor.Magnetometer("Mag",0.01,0,0,0);

%---------------- Environment ------------------------------------------

% Import the environment, the default values give a location of Mojave
% desert with the current date and weather.
env = Env.Environment;

%---------------- Simulation Settings --------------------------------------

% Set the end condition, timestep, simulation fidelity, and outputs
% run 'help Sim.IntegratorSettings' for more details
sim = Sim.IntegratorSettings('full', 0.1, 'medium');

% set the outputs to be shown:
output = 1;

% run rotation visualization (outputs must be on also)
rotationVis = 'on';

% TO-DO: This wind should go into a seperate class script

% change the month for wind data (First 3 letters of month):
month = 'Mar';

% turn wind on and off
windOnOff = 'off';

% create a time array to span the simulation time. Use 500s or more
% w/ recovery on.The code will self-terminate after reaching end
% condition.

if strcmpi('burnout', sim.EndCondition) == 1
    time = constant.burnTime;
elseif ~isnan(str2double(sim.EndCondition)) == 1
    time = round(str2double(endCondition),1);
else
    time = 400;
end

arrayLength = (time / sim.Timestep);
tspan = linspace(0,time,arrayLength+1);

% set the initial position (x,y,z). Accoun ts for starting elevation.
pos = [env.Elevation;0;0];

% set the initial velocity (xdot,ydot,zdot).
vel = [0;0;0];

% initial angle (z angle, y angle, x angle) - following 3-2-1 sequence
angleVector = [0;0.1;0];

% initial rotation rate (x rate, y rate, z rate)
omega = [0;0;0];

% initalize the quaternion based on the euler angle input:
quatVector = eul2quat(angleVector.', "ZYX").';

% initial state vector
Init = [pos;vel;omega;quatVector];

% import aerodynamics data for CMS / R4
if strcmpi(rocket.name, 'CMS') == 1
    rasData = readmatrix("Inputs/RasAeroDataCulled2.CSV");

elseif strcmpi(rocket.name, 'R4') == 1
    rasData = readmatrix("RASAero\Final_with_pumps.CSV");

else
end

% import wind data
%wind = Wind;
%windData = readmatrix("Inputs/WindData.xlsx");
wind = Env.Wind;
windData = wind.windData;
windDataInput = parseWind(windData, month);

% import atmosphere;
atmosphere = readmatrix("Inputs/AtmosphereModel.csv");

% create an array of the center of mass, mass, and moment of inertia of the
% rocket
[totCoM, totMass, MoI] = VariableCoM(sim.Timestep, tspan, 0);

% additional options for RK4 (stop after reaching final condition)
opt = odeset('Events', @(tspan, Init) stoppingCondition(tspan, Init, sim.EndCondition), ...
    'RelTol', sim.relTol, 'AbsTol', sim.absTol);

%---------------- Run the RK4 Integration ----------------------------------
tic;
[timeArray, out] = ode45(@(time,input) RK4Integrator(time,input,rasData,atmosphere, ...
    totCoM,totMass,MoI,windDataInput,windOnOff, rocket, sim), tspan, Init, opt);
toc;

%% Outputs:

% create a struct which contains all of the output information:
outputStruct = struct;

outputStruct.time = timeArray;

% output additional arrays from the integrator
for k = 1:numel(timeArray)
    [~, outputStruct.mach(k,1), outputStruct.AoA(k,1), outputStruct.accel(k,:), ...
        outputStruct.cD(k,:)] = RK4Integrator(timeArray(k), out(k,:), rasData, ...
        atmosphere,totCoM, totMass, MoI, windDataInput, windOnOff, rocket, sim);
end


if output == 1
    % make the outputs real (long monte carlo runs can generate complex values)
    out = real(out);
    outputStruct.AoA = real(outputStruct.AoA);

    % parse rk4 outputs:
    posArray = out(:,1:3);
    velArray = out(:,4:6);
    omega = out(:,7:9)';
    quatArray = out(:,10:13);

    outputStruct.pos = posArray;
    outputStruct.vel = velArray;
    outputStruct.omega = omega;
    outputStruct.quat = quatArray;

    % get the height measurement based on the sensor properties
    heightMeasAltimeter = altimeter.AltitudeMeasurement(posArray(:,1),sim.Timestep, velArray(:,1));

    [posMeasGps, velMeasGps] = gps.GNSSMeasurement(posArray,velArray,sim.Timestep);

    % convert to lat and long for plotting on map:

    E = wgs84Ellipsoid;
    [lats,longs, ~] = ned2geodetic(out(:,3),out(:,2),-out(:,1),env.LatLong(1),env.LatLong(2),E.SemimajorAxis,E);

    % get the outputs from the magnetometer
    xyzMag =  mag.MagnetometerMeasurement(env,[lats,longs,posArray(:,1)], sim.Timestep);

    uif = uifigure;
    g = geoglobe(uif);
    
    geoplot3(g, lats, longs, out(:,1), 'r-', LineWidth= 1)
    campos(g,env.LatLong(1)-0.1,env.LatLong(2)-0.1,15000)
    campitch(g,-30)
    camheading(g,45)
    

    % find end conditions for graphs / animations
    endTime = length(outputStruct.AoA) * sim.Timestep;

    % grab parameters at max Q and off the rail
    [maxVel, maxIndex] = max(out(:,4));
    maxqAccel = outputStruct.accel(maxIndex,1);
    maxqpos = posArray(maxIndex,1);

    machTable = rasData(1:300,1);
    cdTable = rasData(1:300,3);
    maxqMach = outputStruct.mach(maxIndex);
    [~, maxqMachIndex] = min(abs(machTable-maxqMach));
    maxqCD = cdTable(maxqMachIndex);

    [~, railIndex] = min(abs(posArray(1:100,1)-constant.railHeight));
    railVel = out(railIndex,4);
    railAccel = outputStruct.accel(railIndex,1);

    railMach = outputStruct.mach(railIndex);
    [~, railMachIndex] = min(abs(machTable-railMach));
    railCD = cdTable(railMachIndex);

    apogee = max(posArray(:,1));

    fprintf("\nParameters at Max Q:\n")
    fprintf(" Velocity: %.2f m/s\n Mach: %.3f\n Acceleration: %.3f m/s^2\n Drag Coefficient: %.4f\n",maxVel, maxqMach, maxqAccel, maxqCD);
    fprintf("Off-Rail Parameters:\n")
    fprintf(" Velocity: %.2f m/s\n Mach: %.3f\n Acceleration: %.3f m/s^2\n Drag Coefficient: %.4f\n",railVel, railMach, railAccel, railCD);
    fprintf("Rocket Apogee (AMSL): %.2f m\n", apogee)
    fprintf("Rocket Apogee (AGL): %.2f m\n", apogee - env.Elevation)

    

    %% Plotting:

    colorlist = ["#ff595e", "#ff924c", "#ffbe0b", "#8ac926", "#1982c4", "#6a4c93", "#06402B"];

    % Earth Frame XYZ position:
    figure;
    fname = 'Cartesian Elements';

    subplot(2,1,1)
    hold on
    plot(timeArray, posArray(:,1), 'Color', colorlist(1));
    plot(timeArray, posArray(:,2), 'Color', colorlist(2));
    plot(timeArray, posArray(:,3), 'Color', colorlist(3));

    xlim([0, endTime]);
    title("Rocket Position in Earth Frame")
    xlabel("Time (s)")
    ylabel("Position [m]")
    legend("$X$","$Y$","$Z$")
    grid on
    hold off


    subplot(2,1,2)
    hold on
    plot(timeArray, velArray(:,1), 'Color', colorlist(4), 'LineStyle','-');
    plot(timeArray, velArray(:,2), 'Color', colorlist(5), 'LineStyle','-');
    plot(timeArray, velArray(:,3), 'Color', colorlist(6), 'LineStyle','-');
    xlabel("Time (s)")
    title("Rocket Velocity in Earth Frame")
    ylabel("Velocity [m/s]")
    legend("$V_x$", "$V_y$", "$V_z$");
    grid on

    %print(hfig,fname,'-dpdf','-painters','-fillpage')
    %print(hfig,fname,'-dpng','-r00')

    % Euler Angles:
    eulerAngles = quat2eul(quatArray,"ZYX");
    figure;
    plot(timeArray, eulerAngles);
    xlim([0,endTime]);
    title("Euler Angles: 3-2-1")
    xlabel("Time (s)")
    ylabel("Euler Angles")
    legend('psi', 'theta', 'phi');

    % Angle of Attack:
    figure;
    plot(timeArray, outputStruct.AoA);
    xlim([0,endTime]);
    title("Angle of Attack")
    xlabel("Time (s)")
    ylabel("Angle of Attack [deg]")

    % Rocket Trajectory Plot:
    figure;
    plot3(posArray(1:int32(endTime / sim.Timestep),3), posArray(1:int32(endTime / sim.Timestep),2), posArray(1:int32(endTime / sim.Timestep),1))
    % plot3(posArray(1:endTime / dt,3), posArray(1:endTime / dt,2), zeros(endTime / dt), '--')
    % plot3(posArray(1:endTime / dt,3), zeros(endTime / dt), posArray(1:endTime / dt,1), '--')
    % plot3(zeros(endTime / dt), posArray(1:endTime / dt,2), posArray(1:endTime / dt,1), '--')
    view(43,24);
    xlabel('Dist North (m)');
    ylabel('Dist East (m)');
    zlabel('Height (m)');
    axis equal;
    grid minor;

    % Measurements Plot:
    figure;
    plot(timeArray,posArray(:,1),'-')
    hold on
    plot(timeArray,heightMeasAltimeter,'+')
    % fix the errors in this
    plot(timeArray,posMeasGps(:,1), '+')
    xlabel('Time [s]');
    ylabel('Height [m]');
    legend('Simulation', 'Altimeter Measurement', 'Gps Measurement')
    

    if strcmpi('on', rotationVis) == 1
        % run the rotation visualizer script
        playbackSpeed = 3;
        quatArray = quatArray';
        posArray = posArray';

        RotationsVisualizer(posArray, quatArray, timeArray, endTime, sim.Timestep, playbackSpeed, 0);

        %% csv outputs:

        output = horzcat(timeArray, outputStruct.mach);

        writematrix(output, 'Outputs/MachTime.csv')
    end
end