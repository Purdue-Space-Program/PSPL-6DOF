function main_output = Main(rocket, env, settings, rocket_name, make_plots)
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
% non-zero AoA calculations for lift and drag need better modeling
%
% Monte Script is not up to date with main.


%---------------- Sensor Definition ------------------------------------------

% % Make a really bad altimeter for testing
% altimeter = Sensor.Altimeter("Altimeter", 0.25, 20^2,.5, 5, .01);
% 
% % Make a GPS with measurement update:
% gps = Sensor.GNSS("GPS",2, 3^2, .1, 0);
% 
% % Make a magnetometer:
% mag = Sensor.Magnetometer("Mag",0.01,0,0,0);
% 
% % Make a gyroscope:
% gyro = Sensor.Gyroscope("Gyro",.25,1e-4,.005,.01,0);

% create a time array to span the simulation time. Use 500s or more
% w/ recovery on.The code will self-terminate after reaching end
% condition.

if strcmpi('burnout', settings.EndCondition)
    time = rocket.BurnTime;
elseif ~isnan(str2double(settings.EndCondition))
    time = round(str2double(settings.EndCondition),1);
else
    time = 70;
end

arrayLength = (time / settings.Timestep);
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

% import aerodynamics data
rasData = setAeroData(rocket);

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

plot_CoM_graph = 0; % 1 is true, 0 is false
[totCoM, totMass, MoI] = VariableCoM(settings.Timestep, tspan, plot_CoM_graph, rocket, rocket_name);

% additional options for RK4 (stop after reaching final condition)
opt = odeset('Events', @(tspan, Init) stoppingCondition(tspan, Init, settings.EndCondition), ...
    'RelTol', settings.relTol, 'AbsTol', settings.absTol);

%---------------- Run the RK4 Integration ----------------------------------

% this is just for measuring the time it takes to run this, none of the actual simulation happens here - david
% wait i might be wrong - david
% what the fuck is going on here - david
tic;
[timeArray, out] = ode45(@(time,input) RK4Integrator(time,input,atmosphere, ...
    totCoM,totMass,MoI,windData, rocket, settings, env), tspan, Init, opt);
toc;

%% Outputs:

% create a struct which contains all of the output information:
output_struct = struct;
output_struct.time = timeArray;

% output additional arrays from the integrator
for k = 1:numel(timeArray)
    [~, output_struct.mach(k,1), output_struct.AoA(k,1), output_struct.acceleration(k,:), ...
        output_struct.cD(k,:), moment(k,:)] = RK4Integrator(timeArray(k), out(k,:), ...
        atmosphere, totCoM, totMass, MoI, windData, rocket, settings, env);
end


if settings.Outputs == true
    % make the outputs real (long monte carlo runs can generate complex values)
    out = real(out);
    output_struct.AoA = real(output_struct.AoA);

    % parse rk4 outputs:
    position_array = out(:,1:3);
    velocity_array = out(:,4:6);
    omega = out(:,7:9);
    quat_array = out(:,10:13);

    output_struct.position = position_array;
    output_struct.velocity = velocity_array;
    output_struct.omega = omega;
    output_struct.quat = quat_array;
    
    if make_plots == true
        % convert to lat and long for plotting on map:
        E = wgs84Ellipsoid;
        [lats,longs, ~] = ned2geodetic(out(:,3),out(:,2),-out(:,1),env.LatLong(1),env.LatLong(2),E.SemimajorAxis,E);
    
        uif = uifigure;
        g = geoglobe(uif);
        
        geoplot3(g, lats, longs, out(:,1), 'r-', LineWidth= 1)
        campos(g,env.LatLong(1)-0.1,env.LatLong(2)-0.1,15000)
        campitch(g,-30)
        camheading(g,45)
        
        % find end conditions for graphs / animations
        endTime = length(output_struct.AoA) * settings.Timestep;
    end

    % grab parameters at max Q and off the rail
    [max_Q_vertical_velocity, max_Q_index] = max(out(:,4));
    max_Q_index = max_Q_index - 5; % i love euler method!!!!!!!
    max_Q_acceleration = output_struct.acceleration(max_Q_index,1);
    max_Q_altitude = position_array(max_Q_index,1);
    max_Q_horizontal_velocity = norm(output_struct.velocity(max_Q_index, 2:3));

    machTable = rasData(1:300,1);
    cdTable = rasData(1:300,3);
    max_Q_Mach = output_struct.mach(max_Q_index);
    [~, max_Q_Mach_index] = min(abs(machTable-max_Q_Mach));
    max_Q_CD = cdTable(max_Q_Mach_index);

    [~, off_the_rail_index] = min(abs(position_array(1:50,1) - env.railHeight - position_array(1,1)));
    off_the_rail_velocity = out(off_the_rail_index,4);
    off_the_rail_acceleration = output_struct.acceleration(off_the_rail_index,1);

    rail_Mach = output_struct.mach(off_the_rail_index);
    [~, rail_Mach_index] = min(abs(machTable-rail_Mach));
    railCD = cdTable(rail_Mach_index);

    [apogee_altitude, apogee_index] = max(position_array(:,1));

    apogee_horizontal_velocity = norm(output_struct.velocity(apogee_index, 2:3));


    main_output = struct();
    main_output.max_Q_horizontal_velocity = max_Q_horizontal_velocity;
    main_output.max_Q_vertical_velocity = max_Q_vertical_velocity;
    main_output.max_Q_Mach = max_Q_Mach;
    main_output.max_Q_acceleration = max_Q_acceleration;
    main_output.max_Q_CD = max_Q_CD;
    
    main_output.apogee_horizontal_velocity = apogee_horizontal_velocity;
    main_output.apogee_altitude = apogee_altitude - env.Elevation;


    fprintf("\nParameters at Max Q:\n")
    fprintf("\tHorizontal Velocity: %.2f m/s\n", max_Q_horizontal_velocity);
    fprintf("\tAltitude: %.2f m/s\n", max_Q_altitude);
    fprintf("\tVelocity: %.2f m/s\n Mach: %.3f\n Acceleration: %.3f m/s^2\n Drag Coefficient: %.4f\n",max_Q_vertical_velocity, max_Q_Mach, max_Q_acceleration, max_Q_CD);
    fprintf("Off-Rail Parameters:\n")
    fprintf("\tVelocity: %.2f m/s\n Mach: %.4f\n Acceleration: %.3f m/s^2\n Drag Coefficient: %.4f\n",off_the_rail_velocity, rail_Mach, off_the_rail_acceleration, railCD);
    fprintf("\nRocket Apogee (AMSL): %.2f m\n", apogee_altitude)
    fprintf("Rocket Apogee (AGL): %.2f m\n", apogee_altitude - env.Elevation)

    if make_plots == true
    
        %% Plotting:
    
        colorlist = ["#ff595e", "#ff924c", "#ffbe0b", "#8ac926", "#1982c4", "#6a4c93", "#06402B"];
    
        % Earth Frame XYZ position:
        figure;
        fname = 'Cartesian Elements';
    
        subplot(2,1,1)
        hold on
        plot(timeArray, position_array(:,1), 'Color', colorlist(1));
        plot(timeArray, position_array(:,2), 'Color', colorlist(2));
        plot(timeArray, position_array(:,3), 'Color', colorlist(3));
    
        xlim([0, endTime]);
        title("Rocket Position in Earth Frame")
        xlabel("Time (s)")
        ylabel("Position [m]")
        legend("$X$","$Y$","$Z$")
        grid on
        hold off
    
    
        subplot(2,1,2)
        hold on
        plot(timeArray, velocity_array(:,1), 'Color', colorlist(4), 'LineStyle','-');
        plot(timeArray, velocity_array(:,2), 'Color', colorlist(5), 'LineStyle','-');
        plot(timeArray, velocity_array(:,3), 'Color', colorlist(6), 'LineStyle','-');
        xlabel("Time (s)")
        title("Rocket Velocity in Earth Frame")
        ylabel("Velocity [m/s]")
        legend("$V_x$", "$V_y$", "$V_z$");
        grid on
    
        %print(hfig,fname,'-dpdf','-painters','-fillpage')
        %print(hfig,fname,'-dpng','-r00')
    
        % Euler Angles:
        eulerAngles = quat2eul(quat_array,"ZYX");
        figure;
        plot(timeArray, eulerAngles);
        xlim([0,endTime]);
        title("Euler Angles: 3-2-1")
        xlabel("Time (s)")
        ylabel("Euler Angles")
        legend('psi', 'theta', 'phi');
    
    
        % Angle of Attack:
        figure;
        plot(timeArray, output_struct.AoA);
        xlim([0,endTime]);
        title("Angle of Attack")
        xlabel("Time (s)")
        ylabel("Angle of Attack [deg]")
    
        % Rocket Trajectory Plot:
        figure;
        plot3(position_array(1:int32(endTime / settings.Timestep),3), position_array(1:int32(endTime / settings.Timestep),2), position_array(1:int32(endTime / settings.Timestep),1))
        % plot3(posArray(1:endTime / dt,3), posArray(1:endTime / dt,2), zeros(endTime / dt), '--')
        % plot3(posArray(1:endTime / dt,3), zeros(endTime / dt), posArray(1:endTime / dt,1), '--')
        % plot3(zeros(endTime / dt), posArray(1:endTime / dt,2), posArray(1:endTime / dt,1), '--')
        view(43,24);
        xlabel('Dist North (m)');
        ylabel('Dist East (m)');
        zlabel('Height (m)');
        title("Rocket Trajectory")
        axis equal;
        grid minor;
    
        % Moment plot
        figure;
        plot(timeArray, moment)
        xlabel('time (?)');
        ylabel('moments (?)');
        legend('x','y','z')
        title("Rocket Moments")
        
    
        if settings.RotationVis == true
            % run the rotation visualizer script
            playbackSpeed = 3;
            quat_array = quat_array';
            position_array = position_array';
    
            RotationsVisualizer(position_array, quat_array, timeArray, endTime, settings.Timestep, playbackSpeed, 0);
    
            %% csv outputs:
    
            %output = horzcat(timeArray, outputStruct.mach);
    
            %writematrix(output, 'Outputs/MachTime.csv')
        end
    end

    end

end