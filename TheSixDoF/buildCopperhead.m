function rocket = buildCopperhead(cfg, propMasses, aeroPath)
%BUILDCOPPERHEAD  Assemble a Copperhead rocket from a config struct.
%
%  cfg fields (all masses in lbm, lengths in inches):
%    .name
%    .loxLen         LOx tank length
%    .loxShellMass   LOx tank empty shell mass
%    .ethLen         Ethanol tank length
%    .ethShellMass   Ethanol tank empty shell mass
%    .lowerAFLen     Lower airframe length (default 39 in if not set)
%    .midAvionics    Mid avionics mass
%    .raceway        Raceway mass
%    .pumps          Pumps mass (0 for ablative)
%    .chamber        Chamber mass
%    .thrust         Momentum thrust [N]
%    .exitArea       Nozzle exit area [m²]
%    .exitPres       Nozzle exit pressure [Pa]
%    .totalMassFlow  Total propellant mass flow [kg/s]
%
%  propMasses fields (kg):
%    .lox            LOx propellant load
%    .ethanol        Ethanol propellant load

LBM2KG = 0.45359237;
IN2M   = 0.0254;

DIA_M    = 8.84  * IN2M;
WALL_THK = 0.17  * IN2M;

%% ── Fixed section lengths [m] ────────────────────────────────────────────
L.nose    = 36.00 * IN2M;
L.recBay  = 28.50 * IN2M;
L.heBay   = 22.50 * IN2M;
L.upperAF = 23.50 * IN2M;
L.lox     = cfg.loxLen * IN2M;
L.midAF   = 19.50 * IN2M;
L.eth     = cfg.ethLen * IN2M;
if isfield(cfg, 'lowerAFLen') && ~isempty(cfg.lowerAFLen)
    L.lowerAF = cfg.lowerAFLen * IN2M;
else
    L.lowerAF = 39.00 * IN2M;
end
L.inj     =  1.50 * IN2M;
L.chamber = 16.00 * IN2M;

%% ── Positions from nose tip [m] ──────────────────────────────────────────
p.nose    = 0;
p.recBay  = p.nose    + L.nose;
p.heBay   = p.recBay  + L.recBay;
p.upperAF = p.heBay   + L.heBay;
p.lox     = p.upperAF + L.upperAF;
p.midAF   = p.lox     + L.lox;
p.eth     = p.midAF   + L.midAF;
p.lowerAF = p.eth     + L.eth;
p.inj     = p.lowerAF + L.lowerAF;
p.chamber = p.inj     + L.inj;
totalLen  = p.chamber + L.chamber;

%% ── Structural masses [kg] ───────────────────────────────────────────────
M.nose       = 10.00  * LBM2KG;
M.recBay     = 17.00  * LBM2KG;
M.heBay      =  2.00  * LBM2KG;
M.copv       = 15.00  * LBM2KG;
M.upperAF    =  6.00  * LBM2KG;
M.upperPlumb = 20.00  * LBM2KG;
M.loxShell   = cfg.loxShellMass * LBM2KG;
M.midAF      =  7.50  * LBM2KG;
M.midFluids  =  2.00  * LBM2KG;
M.midAvi     = cfg.midAvionics * LBM2KG;
M.ethShell   = cfg.ethShellMass * LBM2KG;
M.raceway    = cfg.raceway * LBM2KG;
M.lowerAF    = 22.50  * LBM2KG;
M.lowerPlumb = 30.00  * LBM2KG;
M.pumps      = cfg.pumps * LBM2KG;
M.inj        =  9.00  * LBM2KG;
M.chamber    = cfg.chamber * LBM2KG;

burnTime = (propMasses.lox + propMasses.ethanol) / cfg.totalMassFlow;

%% ── Rocket object ────────────────────────────────────────────────────────
rocket = Rocket(cfg.name);
rocket.OuterDiameter = DIA_M;
rocket.TotalLength   = totalLen;
rocket.NoseLength    = L.nose;
rocket.NoseGeometry  = 'Von Karman';

dryMass = M.nose + M.recBay + M.heBay + M.copv + M.upperAF + M.upperPlumb + ...
          M.loxShell + M.midAF + M.midFluids + M.midAvi + ...
          M.ethShell + M.raceway + M.lowerAF + M.lowerPlumb + ...
          M.pumps + M.inj + M.chamber;
rocket.TotalMass = dryMass + propMasses.lox + propMasses.ethanol;

%% ── Components ───────────────────────────────────────────────────────────
function addPM(nm, mass, pos, len)
    c = PointMass(nm);
    c.Mass = mass; c.Position = [pos, 0, 0]; c.ComponentLength = len;
    rocket.addComponent(c);
end

% Sections with explicit lengths — distributed over their own span
addPM('Nosecone',      M.nose,       p.nose,    L.nose);
addPM('RecoveryBay',   M.recBay,     p.recBay,  L.recBay);
addPM('HeliumBay',     M.heBay,      p.heBay,   L.heBay);
addPM('UpperAirframe', M.upperAF,    p.upperAF, L.upperAF);
addPM('MidAirframe',   M.midAF,      p.midAF,   L.midAF);
addPM('LowerAirframe', M.lowerAF,    p.lowerAF, L.lowerAF);

% No-length items — distributed over the bay they occupy
addPM('COPV',          M.copv,       p.heBay,   L.heBay);       % dist over Helium Bay
addPM('UpperPlumbing', M.upperPlumb, p.upperAF, L.upperAF);     % dist over Upper Airframe
addPM('MidFluids',     M.midFluids,  p.midAF,   L.midAF);       % dist over Mid Airframe
addPM('MidAvionics',   M.midAvi,     p.midAF,   L.midAF);       % dist over Mid Airframe
addPM('LowerPlumbing', M.lowerPlumb, p.lowerAF, L.lowerAF);     % dist over Lower Airframe
addPM('Raceway',       M.raceway,    p.lox,     L.lox + L.midAF + L.eth);  % dist over LOx + Mid + Eth

if M.pumps > 0
    addPM('Pumps', M.pumps, p.lowerAF, L.lowerAF);              % dist over Lower Airframe
end

t = Tank('LOxTank');
t.Mass = M.loxShell; t.LiquidMass = propMasses.lox;
t.Position = [p.lox,0,0]; t.Length = L.lox;
t.TankDia = DIA_M; t.Thickness = WALL_THK;
rocket.addComponent(t);

t = Tank('EthanolTank');
t.Mass = M.ethShell; t.LiquidMass = propMasses.ethanol;
t.Position = [p.eth,0,0]; t.Length = L.eth;
t.TankDia = DIA_M; t.Thickness = WALL_THK;
rocket.addComponent(t);

e = PropulsionSystem('Engine');
e.Mass            = M.inj + M.chamber;
e.Position        = [p.inj, 0, 0];
e.ComponentLength = L.inj + L.chamber;
e.Thrust          = cfg.thrust;
e.BurnTime        = burnTime;
e.ExitArea        = cfg.exitArea;
e.ExitPressure    = cfg.exitPres;
rocket.addComponent(e);

if ~isempty(aeroPath) && isfile(aeroPath)
    rocket.AeroData = rocket.setAeroData(aeroPath);
end

end
