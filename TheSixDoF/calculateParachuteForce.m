function [drag,paraVector] = calculateParachuteForce(pos, vel, rho, env, tDrogue)
% PSP FLIGHT DYNAMICS:
%
% Title: CalculateParachuteForce
% Author: Preston Wright and Hudson Reynolds - Created: 9/1/25
% Last Modified: 9/9/25
%
% Description: This takes the velocity of the rocket as well as
% environmental parameters and drogue deployment time to model the drag
% from the parachute in a multibody simulation
%
% Inputs: 
% vel - inertial velocity of the rocket [m/s]
% rho - atmospheric density [kg/m^3]
% environment - environment object
% tDrogue - time at which the drogue chute deploys in the simulation [s]
%
% Outputs:
% drag - inertial drag force imparted by the parachute [N]
% paraVector - directional vector along which the parachute is oriented

% ---- Initializations -----------------------------------------------------

% Other Constants
deformVal = 0.70;                       % Deformation value of inflated chute area
                                        % *This is an approximate value due to
                                        % the difficulty to calculate it*

% Parachute parameters
drogueCd = 0.97;                        % cD for the drogue chute
drogueDia = (25/6) * constant.ft2m;     % drogue  diameter [m]

mainCd = 2.2;                           % cD for the main parachute
mainDia = (97/6) * constant.ft2m;       % main chute diameter [m]
mainDeployAlt = 304.8 + env.Elevation;  % main chute deployment altitude [m]

% define paraLength as a function of time:
paraLength = 8*(1-exp(-.5*tDrogue));

%establish the direction of the parachute based on the freestream
%direction:

freestreamDirection = vel / norm(vel); % Normalize the velocity vector to get the direction

% apply the length constraint to the get parachute position vector:
paraVector = -freestreamDirection * paraLength;

% calculate the drag force:

area = deformVal * pi * (0.5 * drogueDia)^2 * (1-exp(-.5*tDrogue));

drag = 1/2 * rho * drogueCd * area * dot(vel,vel);

drag = drag * -freestreamDirection;

end