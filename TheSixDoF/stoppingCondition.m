function [value, isterminal, direction] = stoppingCondition(tspan, Init, condition)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSP FLIGHT DYNAMICS:
%
% Title: VariableCoM
% Author: Preston Wright - Created: 9/28/2024
%
% Description: This function calculates an array for the center of mass in 
% meters from the nose of the rocket. Uses given initialized parameters for 
% the rocket to approximate the center of mass 
%
% Inputs: 
% dt = given simulation time step [s]
% tspan = array of time values for total simulation run time for the given
%         time step [s]
% graph = boolean operator that controls the output of CoM visualizations:
%         1 will output visuals, 0 will not
%
% Outputs:
% totCoM  = array of all center of mass values for the rocket with respect
%           to time [s|m]
% totMass = array of all total mass values for the rocket with respect to
%           time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if strcmpi('apogee', condition) == 1
    value = (Init(4) < 0);
    isterminal = 1;   % Stop the integration
    direction  = 0;
else
    value = (Init(1) < 0);
    isterminal = 1;   % Stop the integration
    direction  = 0;
end