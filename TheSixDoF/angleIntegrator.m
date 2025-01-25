function [eulerRates] = angleIntegrator(omegaVector, angleVector)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSP FLIGHT DYNAMICS:
%
% Title: angleIntegrator
% Authors: Hudson Reynolds and Preston Wright - Created: 9/7/2024
%          Preston Wright - Edited: 9/9/2024
%
% Description: This function takes the omega vector in the body frame and
% the current euler angle and outputs the euler rotation rates
%
% Note: Professor Früh's notes from AAE 340 Section 8.1 and 8.2 go through
% these steps in more detail
%
% Inputs: 
% omegaVector = array of the angular velocity in the body frame [rad/s]
% angleVector = array of the current euler angles [rad]
%
% Outputs:
% eulerRates = the euler rate time derivative in the inertial frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize angles
phi = angleVector(1);       % Roll
theta = angleVector(2);     % Pitch
psi = angleVector(3);       % Yaw

% These are the component rotation matrices
aRollPhi= [1      0            0;
           0  cos(phi)  sin(phi);
           0  -sin(phi) cos(phi)];

aPitchTheta= [cos(theta) 0 -sin(theta);
              0          1            0;
              sin(theta) 0  cos(theta)];

aYawPsi = [cos(psi)  sin(psi) 0;
           -sin(psi) cos(psi) 0;
           0            0     1];

% These are the rotation rate vectors calculated in accordance with our
% rotation order (3-2-1), starting with psi through all frames, then theta
% from intermediate frame 1, and phi as the final rotation to the body 
% frame
psiDot = aRollPhi * aPitchTheta * [0; 0; 1];
thetaDot = aRollPhi * [0; 1; 0];
phiDot = [1; 0; 0];

% This BMatrix is the coefficient matrix of the euler rates calculated in
% accordance to Früh's notes (reference confluence)
bMatrix = [phiDot, thetaDot, psiDot];

eulerRates = inv(bMatrix)*omegaVector;

