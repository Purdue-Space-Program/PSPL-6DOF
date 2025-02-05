%
clear
clc
close all
format longg

pos = [0;0;0];
vel = [0;0;0];
omega = [0;0;0];

Xdata = []; %contains mass
Ydata = []; %contains max height

tic;
for i = -50:50:50
    dt = 0.04;
    maxTime = 50;
    tspan = linspace(0,maxTime,maxTime/dt+1);

    angleVector = [0;0;0];
    quatVector = eul2quat(angleVector.', "XYZ").';
    Init = [pos;vel;omega;quatVector];
    J = eye(3);

    %huge stuff
    rasData = readmatrix("Inputs/Final_with_pumps.CSV");
    windData = readmatrix("Inputs/WindData.xlsx");

    %run the simulation
    mInit = 325 + i;
    [timeArray, out] = ode45(@(time,input) RK4Integrator(time,input,rasData,mInit,J,windData), tspan, Init);

    %collect the results
    posArray = [out(:,1), out(:,2), out(:,3)];

    Xdata(end+1) = mInit;
    Ydata(end+1) = max(posArray(:,1));
    display( [Xdata(end),Ydata(end)] )
end
toc;

plot( Xdata,Ydata, Marker='.',MarkerSize=5,Color='#00aa00' );
xlabel( "Starting dry mass (lbm)", fontSize=18 )
ylabel( "Apogee no-wind (m)", fontSize=18  )
title( "Apogees for different variations in starting mass",fontSize=20 )
grid on