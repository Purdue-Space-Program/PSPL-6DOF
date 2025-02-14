%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSP FLIGHT DYNAMICS:
%
% Title: MainRK4
% Author: Hudson Reynolds
%
% Description: Runs the 6DoF with variable thrust-to-weight-ratio to analyze stability
%
% Inputs: N/A
%
% Outputs:
% see subfunctions for specific outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
clear, clc, close all
TRIALS = 5;

% settings
endCondition = 'full';
windOnOff = 'off'; % <- really slow
rocket = 'CMS';
month = 'Mar';

dt = 0.1;
tspan = 0:dt:500;

%% Initial State
pos = [0;0;0];
vel = [0;0;0];
omega = [0;0;0];
angleVector = [0;0;0];
quatVector = eul2quat(angleVector.', "XYZ").';
Init = [pos;vel;omega;quatVector];

%% Imports
if strcmpi(rocket, 'CMS') == 1
    rasData = readmatrix("Inputs/RasAeroDataCulled.CSV");
elseif strcmpi(rocket, 'R4') == 1
    rasData = readmatrix("RASAero/Final_with_pumps.CSV");
end

windData = readmatrix("Inputs/WindData.xlsx");
windDataInput = parseWind(windData, month);

atmosphere = readmatrix("Inputs/AtmosphereModel.csv");

[totCoM, totMass, MoI] = VariableCoM(dt, tspan, 0);

%output lists
Xdata = [];
Ydata = [];

%% RK4:
opt = odeset( ...
    Events=@(tspan, Init) stoppingCondition(tspan, Init, endCondition), ...
    RelTol=1e-6 ...
);

for trial = 1:TRIALS
    %put whatever needs to be varied in here
    %if the variable is defined in RK4, pass it into the last part 
    %of the function and modify RK4
    ONE_DEG = 0.0174533;
    angleVector = [0;abs(ONE_DEG*20*randn);abs(ONE_DEG*20*randn)];
    quatVector = eul2quat(angleVector.', "XYZ").';
    Init = [pos;vel;omega;quatVector];

    tic;
    [timeArray, out] = ode45(@(time,input) RK4Integrator( ...
        time,input,rasData,atmosphere,totCoM, ...
        totMass,MoI,windDataInput,windOnOff,1), tspan, Init, opt);
    toc;

    ap = max(out(:,1));
    Ydata(end+1)=ap
end

plot(Ydata,marker='.',markersize=25)