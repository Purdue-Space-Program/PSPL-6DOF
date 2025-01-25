%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSP FLIGHT DYNAMICS:
%
% Title: MainRK4
% Author: Hudson Reynolds - Created: 9/21/2024
%
% Description: This is the overarching function that runs the 6-DoF,
% calling all neccesary functions to run the simulation. The overarching
% simulation structure uses an RK4 structure using ODE45
%
% Inputs: N/A
%
% Outputs:
% see subfunctions for specific outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear everything before running the code:
clear
clc
close all

%thrust to weight ratio calcs
g = 9.81;             % gravity constant, in m/s^2.
mInit = 74.69;        % initial mass of the rocket

for i = 2:0.2:8
    TWR = i;
    thrustMag = TWR * mInit * g;
end

% create a time array to span the entire simulation time. Use 500s or more w/ recovery on.
dt = 0.1;
time = 500;
arrayLength = (time / dt);
tspan = linspace(0,time,arrayLength+1);

% position (x,y,z)
pos = [0;0;0];
% velocity (xdot,ydot,zdot)
vel = [0;0;0];    
angleVector = [0;0.1;0.1];
omega = [0;0;0];
%initalize the quaternion based on the euler angle input:
quatVector = eul2quat(angleVector.', "XYZ").';

Init = [pos;vel;omega;quatVector];

%data import
rasData = readmatrix("RasAeroData.CSV");

windData = readmatrix("WindData.xlsx");



%huge CoM and Mass array
[totCoM, totMass] = VariableCoM(dt, tspan, 0);

%run RK4:
tic;
[timeArray, out] = ode45(@(time,input) RK4Integrator(time,input,rasData,totCoM,totMass, windData, 1), tspan, Init);
toc;

for k = 1:numel(timeArray)
    [~, machArray(k,1), AoArray(k,1)] = RK4Integrator(timeArray(k), out(k,:), rasData, totCoM, totMass, windData, 1);
end

% make the outputs real as fuck
out = real(out);
AoArray = real(AoArray);

% Array outputs:
posArray = [out(:,1), out(:,2), out(:,3)];

velArray = [out(:,4), out(:,5), out(:,6)];

omega = [out(:,7), out(:,8), out(:,9)];

quatArray = [out(:,10), out(:,11), out(:,12), out(:,13)];

endTime = min((find(posArray(:,1) < 0, 1)) * dt, arrayLength * dt);

%% Plotting:
% Earth Frame XYZ position:
figure(1)
plot(timeArray, posArray);
xlim([0, endTime]);
title("Rocket XYZ Coordinates in Earth Frame")
xlabel("Time (s)")
ylabel("Distance (m)")
legend("X","Y","Z");

% Earth Frame Rocket Velocity:
figure(2)
plot(timeArray, velArray);
xlim([0, endTime]);
title("Rocket velocity in Earth Frame")
xlabel("Time (s)")
ylabel("Velocity")
legend("X","Y","Z");

% Euler Parameters:
figure(3)
plot(timeArray, quatArray);
xlim([0,endTime]);
title("Euler Parameters")
xlabel("Time (s)")
ylabel("Euler Parameters")
legend('q0', 'q1', 'q2', 'q3');

% Angle of Attack:
figure(4)
plot(timeArray, AoArray);
xlim([0,endTime]);
title("Angle of Attack")
xlabel("Time (s)")
ylabel("Angle of Attack [deg]")

% Rocket Trajectory Plot:

figure(5)
plot3(posArray(1:endTime / dt,3), posArray(1:endTime / dt,2), posArray(1:endTime / dt,1))
% plot3(posArray(1:endTime / dt,3), posArray(1:endTime / dt,2), zeros(endTime / dt), '--')
% plot3(posArray(1:endTime / dt,3), zeros(endTime / dt), posArray(1:endTime / dt,1), '--')
% plot3(zeros(endTime / dt), posArray(1:endTime / dt,2), posArray(1:endTime / dt,1), '--')
view(43,24);
xlabel('Dist North (m)');
ylabel('Dist East (m)');
zlabel('Height (m)');
axis equal;
grid minor;

%% outputs:

output = horzcat(timeArray, machArray);

writematrix(output, 'MachTime.csv')

% run the rotation visualizer script
playbackSpeed = 2;
quatArray = quatArray';
posArray = posArray';

RotationsVisualizer(posArray, quatArray, timeArray, endTime, dt, playbackSpeed, 0);

