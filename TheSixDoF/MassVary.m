%
clear
clc
close all
format longg

pos = [0;0;0];
vel = [0;0;0];
omega = [0;0;0];

Xdata = []; %contains massDry
Ydata = []; %contains max height

tic;
for i = -50:50
    dt = 0.1;
    maxTime = 50;
    tspan = linspace(0,maxTime,maxTime/dt+1);

    angleVector = [0;0;0];
    quatVector = eul2quat(angleVector.', "XYZ").';
    Init = [pos;vel;omega;quatVector];
    J = eye(3);

    %huge stuff
    rasData = readmatrix("Inputs/Final_with_pumps.CSV");

    mInit = 178 + i;
    [totCoM, totMass] = VariableCoM(dt, tspan, mInit, 0);

    windData = readmatrix("Inputs/WindData.xlsx");

    %run the simulation
    [timeArray, out] = ode45(@(time,input) RK4Integrator(time,input,rasData,totCoM,totMass,J,windData), tspan, Init);

    %collect the results
    posArray = [out(:,1), out(:,2), out(:,3)];

    Xdata(end+1) = mInit;
    Ydata(end+1) = max(posArray(:,1));
end
toc;

plot( Xdata,Ydata );