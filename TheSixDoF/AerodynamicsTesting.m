% load the rasAero data

rasAeroData = readmatrix("Inputs/RASAero/RasAeroDataCulled2.csv");

% get the parameters for zero angle of attack

rasMach = rasAeroData(1:300,1);

rasCd = rasAeroData(1:300,3);

rasCl = rasAeroData(1:300,6);

rasCp = rasAeroData(1:300,7); %this is in inches, from the nose

rasCp = rasCp/12 * constant.FT_TO_M; % convert from inches to meters

% load the rocket object and compute its aerodynamics:
rocket = load("Inputs\Saved Rockets\CMS.mat");

rocket = rocket.rocketObj;

% input the rocket into the functions and calculate aerodynamic
% characteristics:

% input a state:
state.alpha = deg2rad(0);
state.rho = 1.225;
state.mu = 1.8e-5;


M = linspace(0,3,100);

for idx = 1:length(M)

    state.Mach = M(idx);
    state.V = M(idx)*340;

    aero(idx) = computeAero(rocket,state);

    drag_model(idx) = aero(idx).CD;
    CP(idx) = aero(idx).x_cp;
    CL(idx) = aero(idx).CN;

end



% plot the drag as a function of mach number:
figure;
plot(rasMach, rasCd);
hold on
plot(M, drag_model)
xlabel('Mach Number');
ylabel('Drag Coefficient ($C_D$)');
title('Drag Coefficient vs Mach Number');
legend('RasAero', 'Ours')
grid on;

% plot the lift as a function of mach number:
figure;
plot(rasMach, rasCl);
hold on
plot(M,CL)
xlabel('Mach Number');
ylabel('Drag Coefficient ($C_L$)');
title('Drag Coefficient vs Mach Number');
grid on;

% plot the center of pressure location as a function of mach number:
figure;
plot(rasMach, rasCp);
hold on
plot(M, CP)
xlabel('Mach Number');
ylabel('Center of Pressure Location (m)');
title('Center of Pressure vs Mach Number');
grid on;



function aero = computeAero(rocket, state)

M = state.Mach;
alpha = state.alpha;

% Precompute useful quantities
%Re = reynoldsNumber(state);
Re = 1e4;
beta = betaFactor(M);

if M < 0.8
    aero = subsonicModel(rocket, state, Re, beta);

elseif M > 1.2
    aero = supersonicModel(rocket, state, beta);

else
    aero_sub = subsonicModel(rocket, state, Re, beta);
    aero_sup = supersonicModel(rocket, state, beta);

    aero = transonicBlend(aero_sub, aero_sup, M);
end

end


function aero = subsonicModel(rocket, state, Re, beta)

alpha = state.alpha;

% --- BODY ---
[CNa_body, x_body] = bodyAero(rocket, 'subsonic');

% --- FINS ---
[CNa_fins, x_fins] = finAero(rocket, beta);

% --- TOTAL NORMAL FORCE ---
CNa = CNa_body + CNa_fins;
CN  = CNa * alpha;

% --- DRAG ---
CD = dragModel(rocket, state, Re);

% --- CENTER OF PRESSURE ---
x_cp = (CNa_body*x_body + CNa_fins*x_fins) / CNa;

% Pack outputs
aero.CN = CN;
aero.CNa = CNa;
aero.CD = CD;
aero.x_cp = x_cp;

end


function aero = transonicBlend(sub, sup, M)

% Smooth sigmoid weight
w = 1 / (1 + exp(20*(M - 1)));

aero.CN  = w*sub.CN  + (1-w)*sup.CN;
aero.CNa = w*sub.CNa + (1-w)*sup.CNa;
aero.CD  = w*sub.CD  + (1-w)*sup.CD;
aero.x_cp = w*sub.x_cp + (1-w)*sup.x_cp;

end

function aero = supersonicModel(rocket, state, beta)

fins = rocket.ComponentList.values{1};

alpha = state.alpha;

% Linear supersonic theory
CNa_body = 4 / beta;
CNa_fins = 4 / beta * fins.Count * 0.5; % simple scaling

x_body = rocket.TotalLength * 0.5;
x_fins = rocket.TotalLength * 0.8;

CNa = CNa_body + CNa_fins;
CN  = CNa * alpha;

CD = dragModel(rocket, state, []);

x_cp = (CNa_body*x_body + CNa_fins*x_fins) / CNa;

aero.CN = CN;
aero.CNa = CNa;
aero.CD = CD;
aero.x_cp = x_cp;

end

function [CNa, x_cp] = bodyAero(rocket, regime)

switch regime
    case 'subsonic'
        CNa = 2;
    otherwise
        CNa = 0;
end

x_cp = 0.5 * rocket.TotalLength;

end


function [CNa, x_cp] = finAero(rocket, beta)

fins = rocket.ComponentList.values{1};

Area = 0.5 * (fins.RootChord + fins.TipChord) * fins.Span;

AR = fins.Span^2 / Area;

CLalpha = (2*pi*AR) / (2 + sqrt(4 + (AR*beta)^2));

CNa = CLalpha * fins.Count;

% Approx CP location
x_cp = fins.Position(1) + fins.RootChord/4;

end

function CD = dragModel(rocket, state, Re)

M = state.Mach;

if isempty(Re)
    Re = 1e6; % fallback
end

Cf = 0.455 / (log10(Re)^2.58);
Cf = Cf / (1 + 0.2*M^2)^0.467;

S_wet = rocket.TotalLength * rocket.OuterDiameter * pi;

Area = pi*(rocket.OuterDiameter/2)^2;

CD_friction = Cf * (S_wet / Area);

% crude base drag
CD_base = 0.03;

% simple wave drag
if M > 1.2
    beta = sqrt(M^2 - 1);
    theta = atan((rocket.OuterDiameter/2)/rocket.NoseLength);
    CD_wave = 2 * theta^2 / beta;
elseif M > 0.8
    CD = 0.3608*M+0.1;
    return

else
    CD_wave = 0;
    CD_friction = CD_friction -1.1;
end

CD = CD_friction + CD_base + CD_wave;

end


function beta = betaFactor(M)
% betaFactor Computes compressibility parameter beta
%
% Input:
%   M     - Mach number
%
% Output:
%   beta  - compressibility factor

    % Small tolerance to avoid sqrt of negative near M = 1
    eps_val = 1e-6;

    if M < 1
        val = 1 - M^2;

        % Clamp to avoid numerical issues near Mach 1
        val = max(val, eps_val);

        beta = sqrt(val);

    elseif M > 1
        val = M^2 - 1;

        val = max(val, eps_val);

        beta = sqrt(val);

    else
        % Exactly Mach 1 (avoid zero)
        beta = sqrt(eps_val);
    end

end