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
    'Inputs/RASAero/r02p02.csv' };

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

%% ── 4 configurations ─────────────────────────────────────────────────────
%
%  cfg(1): r01p02 — Regen,    LOx 32.25in / Eth 39.92in
%  cfg(2): r01p03 — Regen,    LOx 34.14in / Eth 43.67in
%  cfg(3): r02p01 — Ablative, LOx 32.25in / Eth 39.92in
%  cfg(4): r02p02 — Ablative, LOx 34.14in / Eth 43.67in

cfgs(1).name          = 'r01p02';
cfgs(1).loxLen        = 32.25;   cfgs(1).loxShellMass = 13;
cfgs(1).ethLen        = 39.92;   cfgs(1).ethShellMass = 18;
cfgs(1).midAvionics   = 15;
cfgs(1).raceway       = 5;
cfgs(1).pumps         = 16;
cfgs(1).chamber       = 22;
cfgs(1).thrust        = regen.thrust;
cfgs(1).exitArea      = regen.exitArea;
cfgs(1).exitPres      = regen.exitPres;
cfgs(1).totalMassFlow = regen.totalMassFlow;

cfgs(2).name          = 'r01p03';
cfgs(2).loxLen        = 34.14;   cfgs(2).loxShellMass = 15.75;
cfgs(2).ethLen        = 43.67;   cfgs(2).ethShellMass = 21;
cfgs(2).midAvionics   = 15;
cfgs(2).raceway       = 5;
cfgs(2).pumps         = 16;
cfgs(2).chamber       = 22;
cfgs(2).thrust        = regen.thrust;
cfgs(2).exitArea      = regen.exitArea;
cfgs(2).exitPres      = regen.exitPres;
cfgs(2).totalMassFlow = regen.totalMassFlow;

cfgs(3).name          = 'r02p01';
cfgs(3).loxLen        = 32.25;   cfgs(3).loxShellMass = 13;
cfgs(3).ethLen        = 39.92;   cfgs(3).ethShellMass = 18;
cfgs(3).midAvionics   = 5;
cfgs(3).raceway       = 2;
cfgs(3).pumps         = 0;
cfgs(3).chamber       = 5;
cfgs(3).thrust        = abla.thrust;
cfgs(3).exitArea      = abla.exitArea;
cfgs(3).exitPres      = abla.exitPres;
cfgs(3).totalMassFlow = abla.totalMassFlow;

cfgs(4).name          = 'r02p02';
cfgs(4).loxLen        = 34.14;   cfgs(4).loxShellMass = 15.75;
cfgs(4).ethLen        = 43.67;   cfgs(4).ethShellMass = 21;
cfgs(4).midAvionics   = 5;
cfgs(4).raceway       = 2;
cfgs(4).pumps         = 0;
cfgs(4).chamber       = 5;
cfgs(4).thrust        = abla.thrust;
cfgs(4).exitArea      = abla.exitArea;
cfgs(4).exitPres      = abla.exitPres;
cfgs(4).totalMassFlow = abla.totalMassFlow;

%% ── Propellant masses [kg] ───────────────────────────────────────────────
propMasses(1).lox     = 70.80 * LBM2KG;   propMasses(1).ethanol = 64.30 * LBM2KG;   % r01p02 — shorter tanks
propMasses(2).lox     = 84.50 * LBM2KG;   propMasses(2).ethanol = 76.85 * LBM2KG;   % r01p03 — longer tanks
propMasses(3).lox     = 70.80 * LBM2KG;   propMasses(3).ethanol = 64.30 * LBM2KG;   % r02p01 — shorter tanks
propMasses(4).lox     = 84.50 * LBM2KG;   propMasses(4).ethanol = 76.85 * LBM2KG;   % r02p02 — longer tanks

%% ── Build and run ────────────────────────────────────────────────────────
% Set to 1-4: r01p02, r01p03, r02p01, r02p02
RUN = 4;

rocket = buildCopperhead(cfgs(RUN), propMasses(RUN), aeroPaths{RUN});
runSim(rocket, aeroPaths{RUN});
