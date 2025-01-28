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

%
m2ft = 3.28084;

% create a time array to span the entire simulation time. Use 500s or more w/ recovery on.
dt = 0.1;
time = 500;
arrayLength = (time / dt);
tspan = linspace(0,time,arrayLength+1);

% position (x,y,z)
pos = [0;0;0];
% velocity (xdot,ydot,zdot)
svel = 150/m2ft; %starting velocity
vel = [svel*sind(88);svel*cosd(88);0];    
angleVector = [0;0;0];
omega = [0;0;0];
%initalize the quaternion based on the euler angle input:
quatVector = eul2quat(angleVector.', "XYZ").';

Init = [pos;vel;omega;quatVector];

%huge matrix
rasData = readmatrix("arcasData.CSV");

windData = readmatrix("WindData.xlsx");

%huge CoM and Mass array
[totCoM, totMass, MoI] = VariableCoM(dt, tspan, 0);

%run RK4:
tic;
[timeArray, out] = ode45(@(time,input) RK4Integrator(time,input,rasData,totCoM,totMass, MoI, windData, 1), tspan, Init);
toc;

for k = 1:numel(timeArray)
    [~, machArray(k,1), AoArray(k,1), accel(k,:)] = RK4Integrator(timeArray(k), out(k,:), rasData, totCoM, totMass, MoI, windData);
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

colorlist = ["#ff595e", "#ff924c", "#ffbe0b", "#8ac926", "#1982c4", "#6a4c93", "#06402B"];

figure(1)

hfig = figure;  % save the figure handle in a variable

fname = 'Cartesian Elements';

picturewidth = 20; % set this parameter and keep it forever
hw_ratio = .6; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',16) % adjust fontsize to your document
set(hfig,'DefaultLineLineWidth',1)


hold on
plot(timeArray, posArray(:,1), 'Color', colorlist(1));
plot(timeArray, posArray(:,2), 'Color', colorlist(2));
plot(timeArray, posArray(:,3), 'Color', colorlist(3));

xlim([0, 100]);
title("Rocket Cartesian Elements in Earth Frame")
xlabel("Time (s)")
ylabel("Distance [m]")

yyaxis right
plot(timeArray, velArray(:,1), 'Color', colorlist(4), 'LineStyle','-');
plot(timeArray, velArray(:,2), 'Color', colorlist(5), 'LineStyle','-');
plot(timeArray, velArray(:,3), 'Color', colorlist(6), 'LineStyle','-');
ylabel("Velocity[m/s]")
legend("$X$","$Y$","$Z$", "$V_x$", "$V_y$", "$V_z$");

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

grid on

set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(hfig,fname,'-dpdf','-painters','-fillpage')
print(hfig,fname,'-dpng','-r300')

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