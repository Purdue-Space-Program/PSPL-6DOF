function runSim(rocket, aeroPath)
%RUNSIM  Run the 6DoF for a built rocket and write results to Outputs/.
%
%  Behaviour:
%    - Always computes CoM/MoI via VariableCoMMoI if propellant masses and
%      burn time are set (no NaN).
%    - Runs the full 6DoF only when aero data is also available.
%    - Writes a summary text file to Outputs/<rocket.Name>.txt regardless.

hasAero   = ~isempty(aeroPath) && isfile(aeroPath) && ~isempty(rocket.AeroData);
hasMass   = ~isnan(rocket.TotalMass);

% Pull burn time from the Engine component
burnTime = NaN;
comps = values(rocket.ComponentList);
for k = 1:numel(comps)
    if isa(comps{k}, 'PropulsionSystem')
        burnTime = comps{k}.BurnTime;
        break
    end
end
hasBurn = ~isnan(burnTime);

canMassProp = hasMass && hasBurn;
canRunSim   = canMassProp && hasAero;

%% ── Readiness ────────────────────────────────────────────────────────────
fprintf('\n%s\n', rocket.Name)
fprintf('  Wet mass set:   %s\n', yn(hasMass))
fprintf('  Burn time set:  %s\n', yn(hasBurn))
fprintf('  Aero data:      %s\n', yn(hasAero))
fprintf('  → CoM/MoI:      %s\n', yn(canMassProp))
fprintf('  → Full 6DoF:    %s\n\n', yn(canRunSim))

massProps  = [];
simOut     = [];
env        = [];

%% ── Mass properties ──────────────────────────────────────────────────────
if canMassProp
    fprintf('  Computing mass properties...\n')
    [totMass, CoM, MoI, ~] = VariableCoMMoI(rocket);
    massProps = struct('totMass', totMass, 'CoM', CoM, 'MoI', MoI);
end

%% ── Full 6DoF ────────────────────────────────────────────────────────────
if canRunSim
    env = Environment(35.3474, -117.8091);
    env.Elevation = 627.91;
    env.Date      = datetime('now', 'TimeZone', 'UTC');
    env = env.getLocalWeather();
    fprintf('  [Weather] Wind @ 3 km: %.1f m/s  |  Wind @ 10 km: %.1f m/s\n', ...
        norm(interp1(env.Wind(:,1), env.Wind(:,2:3),  3000, 'linear', 'extrap')), ...
        norm(interp1(env.Wind(:,1), env.Wind(:,2:3), 10000, 'linear', 'extrap')));

    s = IntegratorSettings('apogee', 0.1, 'medium');
    s.Outputs     = true;
    s.RotationVis = true;
    s.Wind        = true;

    clear RK4Integrator  % reset persistent drogue state between runs
    fprintf('  Running 6DoF...\n')
    simOut = Main(rocket, env, s);
end

%% ── Write output ─────────────────────────────────────────────────────────
outDir = 'Outputs';
if ~isfolder(outDir); mkdir(outDir); end
outPath = fullfile(outDir, char(rocket.Name) + ".txt");

fid = fopen(outPath, 'w');
fprintf(fid, '%s\n', rocket.Name);
fprintf(fid, 'Generated: %s\n\n', datestr(now, 'yyyy-mm-dd HH:MM'));

fprintf(fid, '── GEOMETRY ─────────────────────────────────\n');
fprintf(fid, '  Total length : %.4f m  (%.2f in)\n', rocket.TotalLength, rocket.TotalLength/0.0254);
fprintf(fid, '  Outer dia    : %.4f m  (%.2f in)\n', rocket.OuterDiameter, rocket.OuterDiameter/0.0254);
fprintf(fid, '  Wet mass     : %.3f kg  (%.2f lbm)\n\n', rocket.TotalMass, rocket.TotalMass/0.45359237);

if ~isempty(massProps)
    fprintf(fid, '── MASS PROPERTIES ──────────────────────────\n');
    fprintf(fid, '  Initial CoM_x (from nose) : %.4f m  (%.2f in)\n', massProps.CoM(1,2), massProps.CoM(1,2)/0.0254);
    fprintf(fid, '  Burnout CoM_x (from nose) : %.4f m  (%.2f in)\n', massProps.CoM(end,2), massProps.CoM(end,2)/0.0254);
    fprintf(fid, '  Initial Ixx               : %.4f kg·m²\n', massProps.MoI(1,2));
    fprintf(fid, '  Initial Iyy               : %.4f kg·m²\n', massProps.MoI(1,3));
    fprintf(fid, '  Initial Izz               : %.4f kg·m²\n', massProps.MoI(1,4));
    fprintf(fid, '  Burn time                 : %.2f s\n\n', burnTime);

    % Write full CoM/MoI time series to CSV
    csvPath = fullfile(outDir, char(rocket.Name) + "_MassProps.csv");
    header = {'t_s','CoM_x_m','CoM_y_m','CoM_z_m','Ixx','Iyy','Izz'};
    data   = [massProps.CoM, massProps.MoI(:,2:4)];
    writematrix(data, csvPath);
    fprintf(fid, '  Full CoM/MoI time series → %s\n\n', csvPath);
else
    fprintf(fid, '[CoM/MoI not computed — fill in propellant masses and burn time]\n\n');
end

if ~isempty(simOut)
    elev = 627.91;

    % Dynamic pressure at each timestep
    speeds   = vecnorm(simOut.vel, 2, 2);
    rho_isa  = interp1(env.Atmosphere(:,1), env.Atmosphere(:,5), simOut.pos(:,3), 'linear', 'extrap');
    q_flight = 0.5 * rho_isa .* speeds.^2;
    [maxQ, maxQIdx] = max(q_flight);
    AoA_at_maxQ = simOut.AoA(maxQIdx);
    t_maxQ      = simOut.time(maxQIdx);

    fprintf(fid, '── 6DoF RESULTS ─────────────────────────────\n');
    apogee_m  = max(simOut.pos(:,3)) - elev;
    fprintf(fid, '  Apogee AGL     : %.1f m  (%.0f ft)\n', apogee_m, apogee_m/0.3048);
    fprintf(fid, '  Max velocity   : %.2f m/s  (%.0f ft/s)\n', max(speeds), max(speeds)/0.3048);
    fprintf(fid, '  Max Mach       : %.3f\n',     max(simOut.mach));
    fprintf(fid, '  Max-Q          : %.1f Pa  (%.2f psi)  (t = %.2f s)\n', maxQ, maxQ/6894.76, t_maxQ);
    fprintf(fid, '  AoA at Max-Q   : %.4f deg\n', AoA_at_maxQ);
else
    fprintf(fid, '[6DoF not run — fill in engine params and provide RASAero CSV]\n');
end

fclose(fid);
fprintf('  Results → %s\n', outPath)

%% ── Save figures ─────────────────────────────────────────────────────────
if ~isempty(simOut) && isfield(simOut, 'figHandles') && ~isempty(simOut.figHandles)
    figDir = fullfile(outDir, char(rocket.Name));
    if ~isfolder(figDir); mkdir(figDir); end
    for k = 1:numel(simOut.figHandles)
        fh = simOut.figHandles{k};
        if ~isvalid(fh); continue; end
        fname = fh.Name;
        if isempty(fname); fname = sprintf('fig%d', k); end
        exportgraphics(fh, fullfile(figDir, [fname '.png']), 'Resolution', 150);
    end
    fprintf('  Figures  → %s/\n', figDir)
end

end

function s = yn(b)
    if b; s = 'YES'; else; s = 'no'; end
end
