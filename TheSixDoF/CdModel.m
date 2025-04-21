function [dragCoeff] = CdModel()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSP FLIGHT DYNAMICS:
%
% Title: CdModel
% Author: Caleb Rice - Created: 2/26/2025
%
% Description: Calculates drag coefficient of the rocket vs mach number using formulas for
%              friction, pressure, and parasitic drag from Barrowman/OpenRocket Thesis
%
% Inputs:
% rocketParams - array, rocket dimensions & parameters, passed from MainRK4
%
% Outputs:
% dragCoeff - array, drag coefficient vs mach #
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clc

%% Inputs

% Rocket Parameters, all units are metric (m, N)
%rocketParams = readmatrix("TheSixDoF\Inputs\CopperheadParameters.xlsx")

bodyDia = 0.219075;
bodyLength = 6.76656;
bodyFineness = bodyLength / bodyDia;

noseconeLength = 1.1;
noseconeType = 2;
noseconeFineness = noseconeLength / bodyDia;

finCount = 3;
finHeight = 0.15;
tipChord = 0.075;
rootChord = 0.5;
finThickness = 0.007;
leadingEdgeProfile = 2;
finShape = 1;

rocketSmoothness = 1; % specifies surface roughness of rocket body

% Environment
kinematicViscAir = 1.5 * (10 ^ -5); % [m2/s]

% Other Parameters
velocity = 0.01:5:680.01;
mach = velocity / 340;

%% Initializations
% Surface Finishes
roughnessValues = [ % [microns]
    0       % perfectly smooth
    0.1     % average glass
    0.5     % polished surfaces
    2       % aircraft-grade sheet metal
    5       % optimum painted surfaces
    15      % planed wooden boards
    20      % mass-production aircraft paint
    50      % smooth cement or steel plating
    100     % asphalt
    150     % dip-galvanized metal
    200     % incorrectly sprayed aircraft paint
    250     % freshly cast iron
    500     % raw wooden boards
    1000    % typical concrete surface
];

noseconeShapes = [
    "Von Karman"
    "Conical"
    "Ogive"
    "Parabolic"
    "Power Series"
];

finProfiles = [
    "Blunt"
    "Rounded"
];

finShapes = [
    "Delta"
    "Trapezoidal - Symmetric"
    "Trapezoidal - Assymmetric"
];

% fin leading edge and trailing edge angles
if finShape < 1 || finShape > 3
    fprintf("ERROR: finShape must be between 1 and 3\n")
elseif finShapes(finShape) == "Delta"
    leadingAngle = atan((rootChord - tipChord) / finHeight);
elseif finShapes(finShape) == "Trapezoidal - Symmetric"
    prompt = "Please enter either the leading or trailing edge angle in radians: ";
    inputAngle = input(prompt);
    trailingAngle = inputAngle;
    leadingAngle = trailingAngle;
elseif finShapes(finShape) == "Trapezoidal - Assymmetric"
    prompt = "Please enter either the leading or trailing edge angle in radians, followed by a space and either LE or TE to indicate which is being specified: ";
    answer = input(prompt, "s");
    [inputAngle, leadingOrTrailing] = strsplit(answer);
    if leadingOrTrailing == "LE"
        leadingAngle = inputAngle;
        leadingEdgeLength = finHeight * tan(leadingAngle);
        trailingAngle = atan((rootChord - leadingEdgeLength - tipChord) / finHeight);
    elseif leadingOrTrailing == "TE"
        trailingAngle = inputAngle;
        trailingEdgeLength = finHeight * tan(trailingAngle);
        leadingAngle = atan((rootChord - trailingEdgeLength - tipChord) / finHeight);
    else
        fprintf("ERROR: the angle measure in radians must be followed by a single space and then either LE or TE and nothing else\n")
    end
end

% fin mean aerodynamic chord (MAC)
finMAC = (2/3) * ((tipChord^3 - rootChord^3) / (tipChord^2 - rootChord^2));

% nose cone lateral area
syms r(d, L, R, K)

if noseconeShapes(noseconeType) == "Von Karman"
    parameter = input("Input a value between 0 and 1/3 for the LV-Haack parameter (0 is VK, 1/3 is haack). k = ", "s");
    r(d, L, R, K) = (R / sqrt(pi))*sqrt(acos(1 - ((2 * d) / L)) - (0.5 * sin(2 * acos(1 - ((2 * d) / L)))) + (K * (sin(acos(1 - (2 * d / L)))) ^ 3));

elseif noseconeShapes(noseconeType) == "Conical"
    parameter = 0;
    r(d, L, R, K) = (d * R) / L;

elseif noseconeShapes(noseconeType) == "Ogive"
    parameter = input("Input a value between 0 and 1 for the ogive parameter (0 is conical, 1 is tangent ogive). k = ", "s");
    r(d, L, R, K) = sqrt((((L^2 + R^2) * (((2 - K) * L)^2 + (K * R)^2)) / (4 * (K * R)^2)) - ((L / K) - d)^2) - sqrt((((L^2 + R^2) * (((2 - K) * L)^2 + (K * R)^2)) / (4 * (K * R)^2)) - (L / K)^2);

elseif noseconeShapes(noseconeType) == "Parabolic"
    parameter = input("Input a value between 0 and 1 for the parabolic parameter (0 is conical, 1 is a full parabola). k = ", "s");
    r(d, L, R, K) = R * (d / L) * ((1 - K * (d / L)) / (2 - K));

elseif noseconeShapes(noseconeType) == "Power Series"
    parameter = input("Input a value between 0 and 1 for the power parameter (0 is a blunt cylinder, 1 is conical). k = ", "s");
    r(d, L, R, K) = R * (d / L)^K;

else
    fprintf("ERROR: invalid nosecone shape, noseconeType must be between 1 and 5\n")
end

r(d) = r(d, noseconeLength, bodyDia/2, parameter);
rPrime(d) = sqrt(1 + (diff(r(d)))^2);
noseconeArea = 2 * pi * vpaintegral(r(d)*rPrime(d), 0, noseconeLength);
noseconeArea = simplify(noseconeArea);

% total lateral area
bodyWetArea = (2 * pi * (bodyDia / 2) * bodyLength) + noseconeArea;
finWetArea = ((rootChord + tipChord) / 2) * finHeight;
wetArea = bodyWetArea + 2 * finCount * finWetArea;

%% Friction Drag
critReynolds = 51 * ((10^-6) * roughnessValues(rocketSmoothness) / bodyLength) ^ -1.039;
reynolds = (velocity .* bodyLength) ./ kinematicViscAir;    % not fully accurate since speed of sound & viscosity will vary with altitude
CF = zeros([1 length(velocity)]);

for x = 1:length(velocity)
    laminarFriction = 1.48 * 10^-2;
    turbulentFriction = (1.5 * log(reynolds(x)) - 5.6) ^ -2;
    roughnessLimitedFriction = 0.032 * (roughnessValues(rocketSmoothness) / bodyLength) ^ 0.2;
    supersonicRoughnessLimited = roughnessLimitedFriction / (1 + 0.18 * mach(x) ^ 2);
    supersonicTurbulent = turbulentFriction / ((1 + 0.15 * mach(x) ^ 2) ^ 0.58);

    if (reynolds(x) <= 10^4) && (reynolds(x) < critReynolds)
        CF(x) = laminarFriction;                                                % laminar friction coeff

    elseif (reynolds(x) > 10^4) && (reynolds(x) <= critReynolds)
        if mach(x) < 1
            CF(x) = turbulentFriction * (1 - 0.1 * mach(x) ^ 2);                % subsonic turbulent friction coeff

        else
            CF(x) = supersonicTurbulent;                                        % supersonic turbulent friction coeff
        end
    elseif (reynolds(x) > critReynolds) && (reynolds(x) > 10^4)
        if mach(x) < 1
            CF(x) = roughnessLimitedFriction * (1 - 0.1 * mach(x) ^ 2);         % subsonic roughness-limited friction coeff

        elseif supersonicTurbulent >= supersonicRoughnessLimited
            CF(x) = supersonicTurbulent;                                        % turbulent used if it is larger than roughness-limited

        else
            CF(x) = supersonicRoughnessLimited;                                 % supersonic roughness-limited friction coeff
        end
    else
        fprintf("ERROR: reynolds number was outside the allowed bounds (not even sure how you did this tbh)\n")
    end
end

CFriction = CF .* (((1 + (1 / (2 * bodyFineness))) * bodyWetArea + (1 + (2 * finThickness / finMAC)) * (2 * finCount * finWetArea)) / wetArea);

%% Pressure Drag
% base drag
baseDrag = zeros([1 length(velocity)]);
stagPress = zeros([1 length(velocity)]);
for x = 1:length(velocity)
    if mach(x) < 1
        baseDrag(x) = 0.12 + 0.13 * mach(x) ^ 2;
        stagPress(x) = 1 + (mach(x) ^ 2 / 4) + (mach(x) ^ 4 / 40);
    else
        baseDrag(x) = 0.25 / mach(x);
        stagPress(x) = 1.84 - (0.76 / mach(x) ^ 2) + (0.166 / mach(x) ^ 4) + (0.035 / mach(x) ^ 6);
    end
end

% fins
CLE = zeros([1 length(mach)]);
CLE_blunt = 0.85 * stagPress;
if leadingEdgeProfile < 1 || leadingEdgeProfile > 2
    fprintf("ERROR: fin type must be either 1 or 2\n")
elseif finProfiles(leadingEdgeProfile) == "Blunt"
    CLE = CLE_blunt;
elseif finProfiles(leadingEdgeProfile) == "Rounded"
    for x = 1:length(mach)
        if mach(x) < 0.9
            CLE(x) = (1 - mach(x) ^ 2) ^ -0.417 - 1;
        elseif mach(x) >= 0.9 & mach(x) < 1
            CLE(x) = 1 - 1.786 * (mach(x) - 0.9);
        else
            CLE(x) = 1.214 - (0.502 / mach(x) ^ 2) + (0.1095 / mach(x) ^ 4) + (0.0231 / mach(x) ^ 6);
        end
    end
end

CLE = CLE * (cos(leadingAngle) ^ 2);

% nose cone
noseconeDragFunc = noseconePressureDrag(bodyDia / 2, noseconeLength, parameter, noseconeType, noseconeFineness, stagPress);
if noseconeDragFunc(1) == 0
    fprintf("ERROR: noseconePressureDrag returned invalid\n")
    quit;
end
noseconePressure = noseconeDragFunc(mach);

% final scaling
finsPressure = (CLE + baseDrag) * (finThickness * finHeight) / (pi * (bodyDia / 2) ^ 2);

CPressure = finsPressure + baseDrag + noseconePressure;

dragCoeff = CPressure + CFriction;

%% Plotting
figure(1)
subplot(2,2,1)
title("Total Drag")
hold on
plot(mach, dragCoeff, "r-")
xlabel("Mach Number")
ylabel("Drag Coefficient")
grid on
hold off

subplot(2,2,2)
title("Nose Cone Pressure Drag")
hold on
plot(mach, noseconePressure, "b-")
xlabel("Mach Number")
ylabel("Drag Coefficient")
grid on
hold off

subplot(2,2,3)
title("Fin Pressure Drag")
hold on
plot(mach, finsPressure, "g-")
xlabel("Mach Number")
ylabel("Drag Coefficient")
grid on
hold off

subplot(2,2,4)
title("Friction vs Pressure Drag")
hold on
plot(mach, CFriction, "k-", mach, CPressure, 'm-')
xlabel("Mach Number")
ylabel("Drag Coefficient")
grid on
hold off

sgtitle("Drag From Various Sources")

saveas(1, 'Drag Graphs.png')