%{
Purdue Space Program - Liquids
StructuresDOF (Determination Of Forces) - main script
Matthew Lewton, edited by Sishir Mahavadi
%}

%% INPUTS 

% VEHICLE DESIGN INPUTS

% Airframe length and mass (discuss with Sishir Mahavadi)
nose = [36, 20]; % [in, lbs]
rec = [28.5, 15];
heTube = [22.5, 16.3]; % [in, lbs]
midAirframe = [19.5, 18]; % [in, lbs]
LoxTnk = [36, 19.6]; % [in, lbs]
mid = [20, 18];
FuelTnk = [34 17.7]; % [in, lbs]  
fincan = [39 60]; % [in, lbs]
engine = [10, 35]; % [in, lbs]

vehicleInputs = [nose; rec; heTube; midAirframe; LoxTnk; mid; FuelTnk; fincan; engine]; % [in, lbs]

% Misc
OD = 8.84; % [in]

% Fins
finNumber = 3; % Number of fins
rootChord = 32; % [in]
tipChord = 7.2; % [in]
sweepLen = 18;
finHeight = 8.6; % [in]
CPnose = convlength(179.1, "in", "m"); % [in] hardcoded from RasAero

% TRAJ INPUTS

V = 1.991e+03; % [ft/s] (velocity at max Q)
velWind = 61;  %[mph] %(Max expected wind gust)
mach = 1.88;   % [mach] (Mach number at max Q)

a = 4.54; % [g] (Acceleration at Max Q)
ag = a;
Thrust = 1800; % [lbf] engine thrust at max q
CD = .46; % Coefficeint of drag at mach 1.67
%rho = 0.8055; % [kg/m^3] Density of air at Burnout @ 19000 ft
[~,~,~,rho] = atmosisa(4794);
%rho = 0.61;
%Q = 1.1795e5; % [Pa] Dynamic pressure at burnout from traj array
Q = 0.5 * rho * (convvel(V, "ft/s", "m/s"))^2;
% Code Inputs
elemLen = .005 ; % [m], length interval for vehicle array

%% SENSORS

SensLoc = sum(vehicleInputs(:,1)) - cumsum(vehicleInputs(:,1)); % [in] (All sensors are at base of part)
SensLabels = ["Nose Cone Bottom" "Rec Bay Bottom" "He Tube Bottom" "Upper Bottom" "Lox Tank Bottom" "Mid Bottom" "Fuel Tank Bottom" "Fin Can Bottom", "Engine Bottom"];

%% Unit Conversion to Metric (fuck slugs)
vehicleInputs(:, 1) = convlength(vehicleInputs(:, 1), "in", "m"); % [m] convert vehicle lengths from in to m
vehicleInputs(:, 2) = convmass(vehicleInputs(:, 2), "lbm", "kg"); % [kg] convert vehicle masses from lbs to kg

rootChord = convlength(rootChord, "in", "m"); % [m] convert to m
tipChord = convlength(tipChord, "in", "m"); % [m] convert to m
sweepLen = convlength(sweepLen, "in", "m"); % [m] convert to m
finHeight = convlength(finHeight, "in", "m"); % [m] convert to m
%finMisalign = finMisalign / 180 * pi; % [rad] convert to rad
%PDistFinMisalign = PDistFinMisalign / 180 * pi; % [rad] convert to rad

V = convvel(V, "ft/s", "m/s"); % [m/s] convert to m/s
velWind = convvel(velWind, "mph", "m/s"); % [m/s] convert to m/s
OD = convlength(OD, "in", "m"); % [m] convert to m

a = convacc(a, "G's", "m/s^2"); %  % [m/s^2] convert to m/s
Thrust = convforce(Thrust, "lbf", "N"); % [N] convert to N

vehicleInputs(:, 1) = elemLen .* (floor(vehicleInputs(:,1) ./ elemLen)); % [m] round all vehicle lengths to the nearest min element size

%% Calculated Parameters

% Vehicle Mass Model
totLen = sum(vehicleInputs(:,1)); % [m]
totMass = sum(vehicleInputs(:,2)); % [kg]

aftZero = []; % [kg] (array of masses at increments)
Lengths = elemLen:elemLen:totLen; % [m] (array of distances from nosecone of corresponding masses in aftZero at even intervals,measures from back of section) 

for x = 1:numel(vehicleInputs(:,1)) %Create array by concatenating arrays for each section
    locDen = elemLen * vehicleInputs(x, 2) / vehicleInputs(x, 1);
    N =  locDen * ones(1, round(vehicleInputs(x, 1) / elemLen));
    aftZero = cat(2, locDen * ones(1, round(vehicleInputs(x, 1) / elemLen)), aftZero); %each array element represents a discrete mass
end

% Misc
S = pi * power(OD/2, 2); % [m^2]
%Q =  (rho * (V^2) / 2) + (59170.90)
%Q = (Thrust - (totMass * a)) / (S * CD); %[pa] Dynamic pressure 

AOA = atan(velWind/V); % [rad] rocket's angle of attack
%AOA = deg2rad(3);
AOAd = rad2deg(AOA); % [deg] rocket's angle of attack
machAdj = 1 / sqrt(mach ^ 2 - 1); %(Divide coefficients calcualted before with Barrowman equations to account for compressible flow)

% Nose
dCLnose =  2 * machAdj; % (Stability derivitive for conical or ogive nose)
cpNose = 0.466 * vehicleInputs(1); % [m] CP of Nose from top of nose cone (assumiung ogive)

% Fins
midChord = sqrt((sweepLen + (tipChord - rootChord) / 2) ^ 2 + (finHeight ^ 2)); % [m] Mid chord of fin
Kfb = (1 + OD / (OD + 2 * finHeight )); % Normal force interference coefficeint
dCLfins = machAdj * Kfb * ( 4 * finNumber * (finHeight / OD) ^2) / (1 + sqrt(1 + (2 * midChord / (rootChord + tipChord)) ^ 2)); % (Stability derivitive for trapezoidal finset)
cpFinsX = ((midChord * (rootChord + 2 * tipChord) / (3 * (rootChord + tipChord)) + (rootChord + tipChord - (rootChord * tipChord / (rootChord + tipChord))) / 6)); % [m] CP of fins from top of root chord

cpFinsY = (finHeight / 3) * (rootChord + 2 * tipChord) / (rootChord + tipChord); % [m] Spanwise CP of fin from root chord
cpFinsYAbs = cpFinsY + OD / 2; % [m] Spanwise CP of fin from longitudinal axis of rocket

finCPAft = rootChord + vehicleInputs(end, 1) - cpFinsX;
finCPNose = totLen - finCPAft;
finArea = finHeight * (rootChord + tipChord) / 2; % [m^2] Area of single fin, used for Torque calculations

% Forces
Lnose = Q * S * AOA * dCLnose; % [N] (Normal force on nose)
Lfins = Q * S * AOA * dCLfins; %[N] (Normal force on fins)

% Center of mass
COMaft = sum(aftZero .* Lengths) / totMass; % [m] from nose 
COMnose =  totLen - COMaft; % [m,] (From aft)

%Center of Pressure
%CPnose = ((dCLnose * cpNose) + ((finCPNose) * dCLfins)) / (dCLfins + dCLnose); % [in] Center of pressure (from nose)

%% AXIAL COMPRESSION
Fc = Thrust - (a * cumsum(aftZero(round(vehicleInputs(end,1) / elemLen):end))); % [N] (axial compression at each inch), truncated top exclude engine
FcLengths = Lengths(round(vehicleInputs(end,1) / elemLen):end); %Truncated lengths array to exclude engine

%% SHEAR
I = sum(aftZero .* ((Lengths - COMaft) .^ 2)); % [kg m^2] rocket's moment of +nertia

ay = (Lnose + Lfins) / totMass; % [m/s^2] (lateral acceleration)
ayn = convacc(ay, "m/s^2", "ft/s^2");
R = (Lfins .* (COMaft - finCPAft) - (Lnose .* (COMnose - cpNose))) / I; % [rad / s^2] (radial acceleration)
Shear = - (ay .* cumsum(aftZero)) - (R * cumsum(aftZero .*  (COMaft - Lengths))); % [N] shear force on each rocket section, Aft is 0!!!!!!

Shear(round((finCPAft) / elemLen):end) = Shear(round((finCPAft) / elemLen):end) + Lfins; %add fin lift above fin cp
Shear(round((totLen-cpNose) / elemLen):end) = Shear(round((totLen-cpNose) / elemLen):end) + Lnose; %nose fin lift above nose cp

%% Bending

Bending = cumsum(Shear) * elemLen; % [N m] Bending moment is the integral of shear,

% %% Torque + Roll Rate
% 
% Tfins = cpFinsYAbs * Q * finMisalign * dCLfins * S; % [N m] Total torque from fin misalignment
% % [rad/s] Terminal roll rate from fin misalignment torque and drag torque balancing out
% rollRate = sqrt((V^2 * finMisalign * dCLfins * S) / (3 * cpFinsYAbs^2 * CDfinNorm * finArea));
% Tdrag = 3/2 * rho * rollRate.^2 * cpFinsYAbs^3 * CDfinNorm * finArea; % [N m] Total torque from counteracting drag normal to fins
% % [rad/s] Terminal roll rate to use for pressure distribution calculations
% PDistRollRate = sqrt((V^2 * PDistFinMisalign * dCLfins * S) / (3 * cpFinsYAbs^2 * CDfinNorm * finArea));
% 
% rollLfin = (Q * PDistFinMisalign * dCLfins * S) / 3; % Lift per fin = Total Lift / 3
% 
% %% Pressure Distributions on Fin from Rolling
% 
% % yAbs = spanwise distance from rocket longitudinal axis
% yAbsStep = linspace(OD / 2, OD / 2 + finHeight, 50);
% 
% % y = spanwise distance from root chord
% yStep = linspace(0, finHeight, 50);
% 
% % Lift distribution using Schrenk approximation
% c = rootChord - (rootChord - tipChord) / finHeight .* yStep;
% ellipseHeight = 4 * finArea / pi / finHeight;
% ellipse = sqrt((1 - yStep.^2 / finHeight^2)) * ellipseHeight;
% 
% Dfin = 1/2 * rho * PDistRollRate.^2 * cpFinsYAbs^2 * CDfinNorm * finArea;

% Check Schrenk approximation
%{
figure();
hold on;
plot(yStep, c);
plot(yStep, ellipse);
hold off;
%}

% [Pa] Lower surface pressure distribution from lift
% PlFin =  1/2 .* (c + ellipse) .* rollLfin / finArea ./ c;

% [Pa] Upper surface (facing direction of roll) pressure distribution from air resistance
% PuFin = 1/2 * rho * PDistRollRate^2 * yAbsStep.^2 * CDfinNorm;

%% Torsional Loads

finAngErr = 0.5; % fin angle error in degrees
finAngErr = finAngErr * pi / 180; % converts from degrees to radians
cL = 2 * pi * finAngErr; % coeff of lift for straight planes
rotForce = cL * Q * (finArea * 3); % lift formula
rotTorq = rotForce * cpFinsYAbs; % torq from the lift force of fins

%% Back to Imperial becuase USA USA USA

totLen = convlength(totLen, "m", "in"); % [in]
totMass = convmass(totMass, "kg", "lbm"); % [lbs]
COMnose = convlength(COMnose, "m", "in"); % [in]
CPnose = convlength(CPnose, "m", "in"); % [in]
OD = convlength(OD, "m", "in"); % [in]
finHeight = convlength(finHeight, "m", "in"); % [in]
% finMisalign = finMisalign / pi * 180; % [deg]
% PDistFinMisalign = PDistFinMisalign / pi * 180; % [deg]

Fc = convforce(Fc, "N", "lbf"); % [lbf] 
Shear = convforce(Shear, "N", "lbf"); % [lbf] 
Lengths = convlength(Lengths, "m", "in"); % [in]
FcLengths = convlength(FcLengths, "m", "in"); % [in]
Bending =  convforce(convlength(Bending, "m", "ft"), "N", "lbf"); % [lbf]
aftZero =  convlength(convforce(aftZero, "N", "lbf"), "m", "in"); % [lb/in]

cpFinsX = convlength(cpFinsX, "m", "in"); % [in]
cpFinsY = convlength(cpFinsY, "m", "in"); % [in]
cpFinsYAbs = convlength(cpFinsYAbs, "m", "in"); % [in]
% yAbsStep = convlength(yAbsStep, "m", "in"); % [in]
% yStep = convlength(yStep, "m", "in"); % [in]
% Tfins = convforce(Tfins, "N", "lbf");
% Tfins = convlength(Tfins, "m", "ft"); % [ft-lb]
% Tdrag = convforce(Tdrag, "N", "lbf");
% Tdrag = convlength(Tdrag, "m", "ft"); % [ft-lb]
% PlFin = convpres(PlFin, "Pa", "psi"); % [psi]
% PuFin = convpres(PuFin, "Pa", "psi"); % [psi]

rotForce = convforce(rotForce, "N", "lbf"); % [lbf]
rotTorq =  convlength(convforce(rotTorq, "N", "lbf"), "m", "ft"); % [lbf/ft]

%% DISPLAY OUTPUTS
fprintf("\nVehicle Parameters:\nTotal Length: %.2f in\nTotal Mass: %.2f lbs\nCOM Location: %.2f in\nCP Location: %.2f in\nStability Margin: %.2f cal\n\nMax Shear: %.2f lbf\nMax Bending Moment: %.2f ft-lbs\nVertical Acceleration: %.2f Gs\nLateral Acceleration: %.2f ft/s^2\nRotational Acceleration: %.2f rad/s^2\nAOA: %.2f deg\n--------------\n", totLen, totMass, COMnose, CPnose, (CPnose-COMnose)/OD, max(abs(Shear)), max(abs(Bending)), ag, ayn, R, AOAd);

% Display forces at sensors
 for i = 1:numel(SensLoc)
    loc = floor(convlength(SensLoc(i), "in", "m") / elemLen) +1 ;
    
    try
        fprintf("%s:\nDistance from aft: %.2f in\nAxial: %.2f lbf\nShear: %.2f lbf\nBending Moment %.2f ft-lbs\n--------------\n", SensLabels(i), Lengths(loc), Fc(floor(convlength(SensLoc(i) - engine(1), "in", "m") / elemLen) +1), Shear(loc), Bending(loc))
    catch
        fprintf("%s:\nDistance from aft: %.2f in\nAxial: %.2f lbf\nShear: %.2f lbf\nBending Moment %.2f ft-lbs\n--------------\n", SensLabels(i), Lengths(loc), -1, Shear(loc), Bending(loc)) %engine bottom is invalid compression output
    end
 end

%fprintf("\nRotational Torq: %f lbf-ft\n", rotTorq)

grapherPSP(1, "Internal Shear Force Along Vehicle Length at Burnout", Lengths, Shear, "Distance from Aft [in]", "Shear Force [lbf]")

grapherPSP(2, "Internal Bending Moment Along Vehicle Length at Burnout", Lengths, Bending, "Distance from Aft [in]", "Bending Moment [ft-lb]")

grapherPSP(3, "Axial Compressive Force Along Vehicle Length at Burnout", FcLengths, Fc, "Distance from Aft [in]", "Compressive Force [lbf]")
 
grapherPSP(4, "Mass Distribution Along Vehicle Length at Burnout", Lengths, aftZero, "Distance from Aft [in]", "Linear Density [lbs/in]")
