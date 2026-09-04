%% Copperhead.m
% Runs all 4 Copperhead configurations.
% Run from repo root:  >> cd TheSixDoF;  Copperhead

clear; clc;
addpath(genpath(fileparts(mfilename('fullpath'))));

% Point MATLAB at the venv that has openmeteo_requests etc.
if ~strcmp(pyenv().Version, '/home/prady/.venvs/sixdof/bin/python3')
    pyenv('Version', '/home/prady/.venvs/sixdof/bin/python3');
end

LBM2KG = 0.45359237;
LBF2N  = 4.44822;
PSI2PA = 6894.76;
Pa_SL  = 101325;

%% ── Aero data ────────────────────────────────────────────────────────────
aeroPaths = { ...
    'Inputs/RASAero/r01p02.csv', ...
    'Inputs/RASAero/r01p03.csv', ...
    'Inputs/RASAero/r02p01.csv', ...
    'Inputs/RASAero/r02p02.csv', ...
    'Inputs/RASAero/r01p01.csv' };

%% ── Regen engine params (Ouroboros) ─────────────────────────────────────
regen.totalMassFlow = (4.30 + 3.91) * LBM2KG;          % kg/s
At_r = regen.totalMassFlow * 1464.66 / (250 * PSI2PA);  % throat area
regen.exitArea  = At_r * 3.81;                           % m²
regen.exitPres  = 11 * PSI2PA;                           % Pa
F_SL_r          = 1500 * LBF2N;                          % N
regen.thrust    = F_SL_r - (regen.exitPres - Pa_SL) * regen.exitArea;

%% ── Ablative engine params (backup) ─────────────────────────────────────
abla.totalMassFlow = (4.065 + 3.695) * LBM2KG;
At_a = abla.totalMassFlow * 1464.48 / (180 * PSI2PA);
abla.exitArea  = At_a * 3.04;
abla.exitPres  = 11 * PSI2PA;
F_SL_a         = 1350 * LBF2N;
abla.thrust    = F_SL_a - (abla.exitPres - Pa_SL) * abla.exitArea;

%% ── 5 configurations ─────────────────────────────────────────────────────
%
%  cfg(1): r01p02 — Regen,    LOx 32.25in / Eth 39.92in / LowerAF 39in
%  cfg(2): r01p03 — Regen,    LOx 38.14in / Eth 47.67in / LowerAF 39in
%  cfg(3): r02p01 — Ablative, LOx 32.25in / Eth 39.92in / LowerAF 39in
%  cfg(4): r02p02 — Ablative, LOx 38.14in / Eth 47.67in / LowerAF 39in
%  cfg(5): r01p01 — Regen,    LOx 36.50in / Eth 32.50in / LowerAF 41in

cfgs(1).name          = 'r01p02';
cfgs(1).loxLen        = 32.25;   cfgs(1).loxShellMass = 22.5;
cfgs(1).ethLen        = 39.92;   cfgs(1).ethShellMass = 24;
cfgs(1).midAvionics   = 15;
cfgs(1).raceway       = 5;
cfgs(1).pumps         = 16;
cfgs(1).chamber       = 22;
cfgs(1).thrust        = regen.thrust;
cfgs(1).exitArea      = regen.exitArea;
cfgs(1).exitPres      = regen.exitPres;
cfgs(1).totalMassFlow = regen.totalMassFlow;

cfgs(2).name          = 'r01p03';
cfgs(2).loxLen        = 38.14;   cfgs(2).loxShellMass = 19;
cfgs(2).ethLen        = 47.67;   cfgs(2).ethShellMass = 21.5;
cfgs(2).midAvionics   = 15;
cfgs(2).raceway       = 5;
cfgs(2).pumps         = 16;
cfgs(2).chamber       = 22;
cfgs(2).thrust        = regen.thrust;
cfgs(2).exitArea      = regen.exitArea;
cfgs(2).exitPres      = regen.exitPres;
cfgs(2).totalMassFlow = regen.totalMassFlow;

cfgs(3).name          = 'r02p01';
cfgs(3).loxLen        = 32.25;   cfgs(3).loxShellMass = 22.5;
cfgs(3).ethLen        = 39.92;   cfgs(3).ethShellMass = 24;
cfgs(3).midAvionics   = 5;
cfgs(3).raceway       = 2;
cfgs(3).pumps         = 0;
cfgs(3).chamber       = 5;
cfgs(3).thrust        = abla.thrust;
cfgs(3).exitArea      = abla.exitArea;
cfgs(3).exitPres      = abla.exitPres;
cfgs(3).totalMassFlow = abla.totalMassFlow;

cfgs(4).name          = 'r02p02';
cfgs(4).loxLen        = 38.14;   cfgs(4).loxShellMass = 19;
cfgs(4).ethLen        = 47.67;   cfgs(4).ethShellMass = 21.5;
cfgs(4).midAvionics   = 5;
cfgs(4).raceway       = 2;
cfgs(4).pumps         = 0;
cfgs(4).chamber       = 5;
cfgs(4).thrust        = abla.thrust;
cfgs(4).exitArea      = abla.exitArea;
cfgs(4).exitPres      = abla.exitPres;
cfgs(4).totalMassFlow = abla.totalMassFlow;

cfgs(5).name          = 'r01p01';
cfgs(5).loxLen        = 36.50;   cfgs(5).loxShellMass = 16.5;
cfgs(5).ethLen        = 32.50;   cfgs(5).ethShellMass = 16;
cfgs(5).lowerAFLen    = 41.00;
cfgs(5).midAvionics   = 15;
cfgs(5).raceway       = 5;
cfgs(5).pumps         = 16;
cfgs(5).chamber       = 22;
cfgs(5).thrust        = regen.thrust;
cfgs(5).exitArea      = regen.exitArea;
cfgs(5).exitPres      = regen.exitPres;
cfgs(5).totalMassFlow = regen.totalMassFlow;

%% ── Propellant masses [kg] ───────────────────────────────────────────────
propMasses(1).lox     = 70.80 * LBM2KG;   propMasses(1).ethanol = 64.30 * LBM2KG;   % r01p02
propMasses(2).lox     = 84.50 * LBM2KG;   propMasses(2).ethanol = 76.85 * LBM2KG;   % r01p03
propMasses(3).lox     = 70.80 * LBM2KG;   propMasses(3).ethanol = 64.30 * LBM2KG;   % r02p01
propMasses(4).lox     = 84.50 * LBM2KG;   propMasses(4).ethanol = 76.85 * LBM2KG;   % r02p02
propMasses(5).lox     = 84.50 * LBM2KG;   propMasses(5).ethanol = 64.30 * LBM2KG;   % r01p01

%% ── Build and run ────────────────────────────────────────────────────────
% Set RUN to a scalar to run one config, or a vector to run multiple.
% Config 5 (r01p01) is the baseline — excluded by default.
RUN = 1:4;
for r = RUN
    rocket = buildCopperhead(cfgs(r), propMasses(r), aeroPaths{r});
    runSim(rocket, aeroPaths{r});
end
