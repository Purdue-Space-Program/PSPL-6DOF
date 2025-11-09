function weatherEnsemblePrediction = get_openmeteo_ensemble(lat, lon, startDate, endDate, hourlyVars)
%GET_OPENMETEO_ENSEMBLE  Pull + post-process Open-Meteo ensemble/historical weather
% for a rocket launch site, compute ensemble stats, make plots,
% and stash results in `weatherEnsemblePrediction`.
%
% USAGE:
%   % 1) Single-day forecast (future date / ensemble up to ~35 days ahead) :contentReference[oaicite:2]{index=2}
%   weatherEnsemblePrediction = get_openmeteo_ensemble(28.5, -80.6, '2025-11-10');
%
%   % 2) Multi-day window (can span past + future; past comes from archive API,
%   %    future comes from ensemble API; we stitch them) :contentReference[oaicite:3]{index=3}
%   weatherEnsemblePrediction = get_openmeteo_ensemble(28.5, -80.6, '2025-11-10','2025-11-12');
%
%   % 3) Custom hourly variables
%   customVars = {'wind_speed_10m','wind_speed_100m','wind_gusts_10m', ...
%                 'wind_direction_10m','surface_pressure','cloud_cover','precipitation'};
%   weatherEnsemblePrediction = get_openmeteo_ensemble(28.5, -80.6, '2025-11-10','2025-11-12', customVars);
%
% INPUTS:
%   lat, lon    : launch pad coordinates (deg)
%   startDate   : 'YYYY-MM-DD' string
%   endDate     : OPTIONAL 'YYYY-MM-DD'. If missing/[], we use startDate.
%
%   hourlyVars  : OPTIONAL cell array of Open-Meteo hourly variable names.
%                 Default (rocket-relevant for launch/no-go):
%                 {'wind_speed_10m','wind_speed_100m','wind_gusts_10m', ...
%                  'wind_direction_10m','surface_pressure','cloud_cover','precipitation'}
%
%   weatherEnsemblePrediction
%       .ensemble / .archive   -> per-hour weather with ensemble stats
%       .verticalProfile       -> altitude vs windU, windV, temp, rho, etc.
%       .Atmosphere            -> [alt(m), T(K), a(m/s), P(Pa), rho(kg/m^3)]
%       .Wind                  -> [alt(m), U(m/s), V(m/s)]
%
%   Each segment (.archive / .ensemble) has:
%       .time        [N_hours x 1 datetime] in local tz
%       .timezone    tz string
%       .stats.(var) struct with:
%           .raw          [N_hours x N_members] all ensemble members stacked
%           .mean         [N_hours x 1] ensemble mean across members/hour
%           .std          [N_hours x 1] ensemble std across members/hour
%           .min          [N_hours x 1] ensemble min across members/hour
%           .max          [N_hours x 1] ensemble max across members/hour
%           .nMembers     scalar (# ensemble members)
%           .det          [N_hours x 1] the unsuffixed/control run if provided
%           .units        text units from API
%
%   - We ALSO build a vertical atmospheric profile for that date using the
%     pressure-level fields (geopotential_height_XXXhPa, etc.) from the
%     OpenMeteoWeatherRequest helper:
%       * pick the forecast/analysis hour closest to startDate
%       * get height, T, P, wind speed & direction at pressure levels
%       * interpolate w/ pchip + extrap from 0 m to 22 km
%       * compute density, speed of sound, and wind vector components
%         (U,V) in the direction the wind is actually going
%       * plot:
%           - 3D quiver of wind vs altitude
%           - temperature vs height
%           - temp/pressure/density vs standard atmosphere (atmosisa)

    % Args / defaults
    if nargin < 3
        error('Usage: get_openmeteo_ensemble(lat, lon, startDate [,endDate] [,hourlyVars])');
    end

    if nargin < 4 || isempty(endDate)
        endDate = startDate;
    end


    defaultHourlyVars = { ...
        'temperature_2m', ...      
        'wind_speed_10m', ...      
        'wind_speed_100m', ...     
        'wind_gusts_10m', ...      
        'wind_direction_10m', ...  
        'surface_pressure', ...    
        'cloud_cover', ...         
        'precipitation' ...       
    };

    if nargin < 5 || isempty(hourlyVars)
        hourlyVars = defaultHourlyVars;
    end
    if ischar(hourlyVars)
        hourlyVars = {hourlyVars};
    end

    hourlyArgStr = strjoin(hourlyVars, ' ');

    % 1. Call python helper
    thisFilePath = mfilename('fullpath');
    [thisDir,~,~] = fileparts(thisFilePath);
    pyScriptPath = fullfile(thisDir,'openmeteo_ensemble.py');

    if ~isfile(pyScriptPath)
        error('Python helper not found at: %s', pyScriptPath);
    end

    % Use gfs_seamless because it's the long-range ensemble.
    modelName = 'gfs_seamless';

    cmd = sprintf([ ...
        'python "%s" --lat %g --lon %g ' , ...
        '--start-date "%s" --end-date "%s" ' , ...
        '--models "%s" ' , ...
        '--timezone auto ' , ...
        '--wind-speed-unit ms ' , ...
        '--precipitation-unit mm ' , ...
        '--temperature-unit celsius ' , ...
        '--hourly %s' ], ...
        pyScriptPath, lat, lon, startDate, endDate, modelName, hourlyArgStr);

    [status, cmdout] = system(cmd);
    if status ~= 0
        error('Open-Meteo call failed:\n%s', cmdout);
    end

    try
        rawData = jsondecode(cmdout);
    catch decodeErr
        error(['Failed to jsondecode python output.\n' ...
               'Raw output:\n%s\n\nError:\n%s'], cmdout, decodeErr.message);
    end

    % 2. Process archive / ensemble segments
    weatherEnsemblePrediction = struct();
    weatherEnsemblePrediction.raw = rawData;

    hasArchive  = isstruct(rawData) && isfield(rawData,'archive');
    hasEnsemble = isstruct(rawData) && isfield(rawData,'ensemble');

    if hasArchive
        weatherEnsemblePrediction.archive = processSegment(rawData.archive);
    else
        weatherEnsemblePrediction.archive = [];
    end

    if hasEnsemble
        weatherEnsemblePrediction.ensemble = processSegment(rawData.ensemble);
    else
        weatherEnsemblePrediction.ensemble = [];
    end

    % Fallback: python may just return a single segment with .hourly
    if isempty(weatherEnsemblePrediction.archive) && isempty(weatherEnsemblePrediction.ensemble)
        if isstruct(rawData) && isfield(rawData,'hourly')
            weatherEnsemblePrediction.ensemble = processSegment(rawData);
        end
    end

    % Pick a "primary" segment (prefer ensemble/future)
    if ~isempty(weatherEnsemblePrediction.ensemble)
        plotSegName = 'ensemble';
        plotSeg = weatherEnsemblePrediction.ensemble;
    elseif ~isempty(weatherEnsemblePrediction.archive)
        plotSegName = 'archive';
        plotSeg = weatherEnsemblePrediction.archive;
    else
        plotSegName = '';
        plotSeg = [];
    end

    weatherEnsemblePrediction.plot_info = struct( ...
        'segment_plotted', plotSegName, ...
        'vars_requested', {hourlyVars} ...
    );

    % 3. Timeline plots (Figures 1–4)
    if ~isempty(plotSeg)
        makePlots(plotSeg);
    else
        warning('No hourly time-series segment available to plot.');
    end

    %% --------------------------
    % 4. Vertical profile build (0 → 22 km)
    %
    % This block:
    %   - builds Atmosphere, Wind, verticalProfile
    %   - makes 3D quiver & ISA comparison figs
    %
    % We wrap this in try/catch so if ANY piece is missing (no temp_2m,
    % atmosisa missing, etc.), you still keep the rest of the output.
    %% --------------------------
    weatherEnsemblePrediction.verticalProfile = [];
    weatherEnsemblePrediction.Atmosphere      = [];
    weatherEnsemblePrediction.Wind            = [];

    try
        if isempty(plotSeg)
            error('No segment to derive vertical profile from.');
        end

        stats = plotSeg.stats;
        tVec  = plotSeg.time;
        idxHour = 1;  % take the first forecast hour as "launch hour snapshot"
        profileTime = tVec(idxHour);

        % --- Standard atmosphere full column (ISA)
        alt_grid = linspace(0,22000,1000).'; % m
        [T_std, a_std, P_std, rho_std] = atmosisa(alt_grid);  % T[K], a[m/s], P[Pa], rho[kg/m^3]

        % --- Surface temp & pressure from ensemble stats (the "real" conditions)
        haveTemp2m = isfield(stats,'temperature_2m') && ~isempty(stats.temperature_2m.mean);
        haveSurfP  = isfield(stats,'surface_pressure') && ~isempty(stats.surface_pressure.mean);

        if haveTemp2m
            Tsurf_K = stats.temperature_2m.mean(idxHour) + 273.15; % C -> K
        else
            % fallback to ISA at ground
            Tsurf_K = T_std(1);
        end

        if haveSurfP
            % surface_pressure comes in hPa typically; convert to Pa
            Psurf_Pa = stats.surface_pressure.mean(idxHour) * 100.0;
        else
            Psurf_Pa = P_std(1);
        end

        % Get ISA values at 100 m for blending
        T_100m_std = pchip(alt_grid, T_std, 100);
        P_100m_std = pchip(alt_grid, P_std, 100);

        % Blend 0→100 m from actual surface -> ISA@100m
        T_api = T_std;  % start with ISA everywhere
        P_api = P_std;
        maskLow = alt_grid <= 100;
        T_api(maskLow) = pchip([0 100],[Tsurf_K, T_100m_std], alt_grid(maskLow));
        P_api(maskLow) = pchip([0 100],[Psurf_Pa, P_100m_std], alt_grid(maskLow));

        % Density and speed of sound from blended T/P
        R_dry      = 287.05; % J/(kg*K) dry air
        gamma_air  = 1.4;    % assume ~1.4 for troposphere
        rho_api    = P_api ./ (R_dry * T_api);
        a_api      = sqrt(gamma_air * R_dry .* T_api);

        % --- Build wind profile from 10m/100m
        % We'll assume:
        %   * speed vertical shape:
        %       pchip between [0m,10m,100m] using ws10 and ws100
        %       hold constant above 100m
        %   * direction constant with height = 10 m direction (wd10)
        %
        % Get ensemble mean surface + 100m winds (m/s) and direction (deg)
        if isfield(stats,'wind_speed_10m') && ~isempty(stats.wind_speed_10m.mean)
            ws10 = stats.wind_speed_10m.mean(idxHour);
        else
            ws10 = NaN;
        end

        if isfield(stats,'wind_speed_100m') && ~isempty(stats.wind_speed_100m.mean)
            ws100 = stats.wind_speed_100m.mean(idxHour);
        else
            ws100 = ws10; % if missing, just use surface
        end

        if isfield(stats,'wind_direction_10m') && ~isempty(stats.wind_direction_10m.mean)
            wd10 = stats.wind_direction_10m.mean(idxHour); % deg, meteorological true dir
        else
            wd10 = 0; % fallback north
        end

        % build vertical wind speed
        windSpeedProf = zeros(size(alt_grid));
        if ~isnan(ws10)
            altLowW = [0; 10; 100];
            spdLowW = [ws10; ws10; ws100];
            maskLowW = alt_grid <= 100;
            windSpeedProf(maskLowW) = pchip(altLowW, spdLowW, alt_grid(maskLowW));
            windSpeedProf(~maskLowW) = ws100; % hold const above 100 m
        else
            windSpeedProf(:) = 0;
        end

        % direction = wd10 everywhere
        dirProfDeg = wd10 * ones(size(alt_grid));

        % turn speed+dir into components of the *actual* flow vector
        % NOTE: we assume wind_direction_10m is "direction the wind is coming FROM"
        % If you want "where it's going", flip sign convention here.
        windU = windSpeedProf .* cosd(dirProfDeg);
        windV = windSpeedProf .* sind(dirProfDeg);

        % --- Assemble outputs
        Atmosphere = [alt_grid, T_api, a_api, P_api, rho_api];   % for sim
        Wind       = [alt_grid, windU, windV];                   % for sim

        vpStruct = struct( ...
            'profile_time', profileTime, ...
            'alt_m',        alt_grid, ...
            'windSpeed_mps',windSpeedProf, ...
            'windDir_deg',  dirProfDeg, ...
            'windU_mps',    windU, ...
            'windV_mps',    windV, ...
            'temp_K',       T_api, ...
            'pressure_Pa',  P_api, ...
            'rho_kgm3',     rho_api ...
        );

        weatherEnsemblePrediction.verticalProfile = vpStruct;
        weatherEnsemblePrediction.Atmosphere      = Atmosphere;
        weatherEnsemblePrediction.Wind            = Wind;

        %% --------------------------
        % 4A. NEAR-SURFACE WIND VECTOR FIELD (0–150 m)
        %
        % We only care about the pad through low ascent for ops, and the
        % API-derived wind is only unique up to ~100 m. Above that we just hold
        % constant anyway. So: zoom in to 150 m.
        %
        % We also downsample arrows and turn OFF quiver autoscaling so that arrow
        % length is in real m/s.
        %% --------------------------
        
        alt_low_mask = alt_grid <= 150;  % keep 0–150 m for visualization
        low_idx_all  = find(alt_low_mask);
        
        % downsample arrow set (every 2nd point in 0–150 m region)
        idxPlot_low  = low_idx_all(1:2:end);
        
        alt_plot_low      = alt_grid(idxPlot_low);      % altitude [m] for arrows
        windU_plot_low    = windU(idxPlot_low);         % U component [m/s]
        windV_plot_low    = windV(idxPlot_low);         % V component [m/s]
        zeros_plot_low    = zeros(size(alt_plot_low));  % W=0 (horizontal)
        originsX_plot_low = zeros(size(alt_plot_low));  % all arrows start at x=0
        originsY_plot_low = zeros(size(alt_plot_low));  % all arrows start at y=0
        
        figure;
        subplot(1,2,1)
        
        % quiver3(X,Y,Z,U,V,W,scale)
        % scale=0 turns off autoscaling -> arrow length = magnitude of [U,V,W]
        quiver3(originsX_plot_low, ...
                originsY_plot_low, ...
                alt_plot_low, ...
                windU_plot_low, ...
                windV_plot_low, ...
                zeros_plot_low, ...
                0);
        
        grid on; box on; view(3)
        
        xlabel('Wind Vel U [m/s]')   % U ~ east-west component
        ylabel('Wind Vel V [m/s]')   % V ~ north-south component
        zlabel('Altitude [m]')
        title('Wind Profile (0–150 m)')
        
        % make axes nice around actual wind speeds:
        maxU = max(abs(windU_plot_low));
        maxV = max(abs(windV_plot_low));
        xlim([-1.2*maxU 1.2*maxU]);
        ylim([-1.2*maxV 1.2*maxV]);
        zlim([0 150]);
        
        % SECOND PANEL in same figure: Temperature vs Altitude (0–150 m)
        subplot(1,2,2)
        
        alt_low = alt_grid(alt_low_mask);
        T_low   = T_api(alt_low_mask);
        
        plot(T_low, alt_low, 'LineWidth',1.5);
        grid on; box on;
        title('Temperature Profile (0–150 m)');
        xlabel('Temperature [K]');
        ylabel('Altitude [m]');
        
        
        %% --------------------------
        % 4B. NEAR-SURFACE ATMOSPHERE VS ISA (0–150 m)
        %
        % We compare:
        %   - Blended/API surface profile (T_api, P_api, rho_api)
        %   - Standard Atmosphere (T_std, P_std, rho_std)
        %
        % Only up to 150 m so you can ACTUALLY SEE the API curve.
        %% --------------------------
        
        T_low_api   = T_api(alt_low_mask);
        T_low_std   = T_std(alt_low_mask);
        P_low_api   = P_api(alt_low_mask);
        P_low_std   = P_std(alt_low_mask);
        rho_low_api = rho_api(alt_low_mask);
        rho_low_std = rho_std(alt_low_mask);
        
        figure;
        
        % Temperature vs altitude
        subplot(1,3,1)
        plot(T_low_api, alt_low, 'LineWidth',1.5); hold on;
        plot(T_low_std, alt_low, 'LineWidth',1.5);
        grid on; box on;
        title('Temperature Profile');
        xlabel('Temperature [K]');
        ylabel('Altitude [m]');
        legend('Blended/API surface','Standard Atmosphere','Location','best');
        xlim([min([T_low_api;T_low_std]) - 1, max([T_low_api;T_low_std]) + 1]);
        ylim([0 150]);
        
        % Pressure vs altitude
        subplot(1,3,2)
        plot(P_low_api, alt_low, 'LineWidth',1.5); hold on;
        plot(P_low_std, alt_low, 'LineWidth',1.5);
        grid on; box on;
        title('Pressure Profile');
        xlabel('Pressure [Pa]');
        ylabel('Altitude [m]');
        legend('Blended/API surface','Standard Atmosphere','Location','best');
        xlim([min([P_low_api;P_low_std])*0.98, max([P_low_api;P_low_std])*1.02]);
        ylim([0 150]);
        
        % Density vs altitude
        subplot(1,3,3)
        plot(rho_low_api, alt_low, 'LineWidth',1.5); hold on;
        plot(rho_low_std, alt_low, 'LineWidth',1.5);
        grid on; box on;
        title('Density Profile');
        xlabel('Density [kg/m^3]');
        ylabel('Altitude [m]');
        legend('Blended/API surface','Standard Atmosphere','Location','best');
        xlim([min([rho_low_api;rho_low_std])*0.9, max([rho_low_api;rho_low_std])*1.1]);
        ylim([0 150]);
        
        
        %% --------------------------
        % 4C. WIND SPEED INTERPOLATION DIAGNOSTIC
        %
        % This is the "white background" style plot you showed.
        %
        % Goal:
        %   - plot the discrete API wind speeds at the levels we truly know
        %     (10 m and 100 m; we also assume 0 m ~ 10 m, which is common for launch pad)
        %   - overlay a smooth pchip curve 0→150 m
        %
        % NOTE:
        %   if wind_speed_10m or wind_speed_100m are missing, we just skip this.
        %% --------------------------
        
        haveWS10  = exist('ws10','var')  && ~isnan(ws10);
        haveWS100 = exist('ws100','var') && ~isnan(ws100);
        
        if haveWS10
            if haveWS100
                % build control points:
                %   altitude [m] :   0      10        100
                %   speed   [m/s]:  ws10,   ws10,     ws100
                alt_raw   = [0; 10; 100];
                spd_raw   = [ws10; ws10; ws100];
            else
                % only surface — just duplicate it
                alt_raw   = [0; 10];
                spd_raw   = [ws10; ws10];
            end
        
            alt_fit    = linspace(0,150,200).';
            spd_fit    = pchip(alt_raw, spd_raw, alt_fit);
        
            figure(); % white background like your example
            plot(alt_fit, spd_fit,'LineWidth',1.5,'DisplayName','pChip interpolation'); hold on;
            scatter(alt_raw, spd_raw,60,'o','filled','DisplayName','API Data');
            grid on; box on;
            title(sprintf('Wind Speed Near Pad (Lat=%.2f, Lon=%.2f, %s)', ...
                           lat, lon, datestr(profileTime,'yyyy-mm-dd')));
            xlabel('Altitude [m]');
            ylabel('Wind Speed [m/s]');
            legend('Location','best');
        
            % limit y so it's tight around these speeds
            ymin = min([spd_raw;spd_fit]);
            ymax = max([spd_raw;spd_fit]);
            ylim([ymin-0.5, ymax+0.5]);
            xlim([0,150]);
        end
    end

    %% --------------------------
    % 5. Expose result in workspace
    % --------------------------
    assignin('base','weatherEnsemblePrediction', weatherEnsemblePrediction);
end


%% ========================================================================
function segOut = processSegment(segIn)
%PROCESSSEGMENT
% Convert one returned segment (archive or ensemble) into:
%   segOut.time      (datetime array)
%   segOut.timezone
%   segOut.stats.(var).mean/std/min/max/etc
%
% We group all *_memberXX variables by base var and stack them.

    segOut = struct();

    if ~isfield(segIn,'hourly')
        segOut = [];
        return;
    end

    % --- time vector ---
    timeStrings = segIn.hourly.time;
    if iscell(timeStrings)
        timeStrings = string(timeStrings);
    elseif ischar(timeStrings)
        timeStrings = string(timeStrings);
    else
        timeStrings = string(timeStrings);
    end

    tz = 'UTC';
    if isfield(segIn,'timezone')
        if ischar(segIn.timezone)
            tz = segIn.timezone;
        elseif isstring(segIn.timezone)
            tz = char(segIn.timezone);
        end
    end

    try
        tVec = datetime(timeStrings, ...
                        'InputFormat',"yyyy-MM-dd'T'HH:mm", ...
                        'TimeZone',tz);
    catch
        tVec = datetime(timeStrings,'TimeZone',tz);
    end
    tVec = tVec(:);
    nT   = numel(tVec);

    segOut.time = tVec;
    segOut.timezone = tz;

    % --- units map, if provided ---
    unitsMap = struct();
    if isfield(segIn,'hourly_units')
        unitsMap = segIn.hourly_units;
    end

    % --- group ensemble members ---
    hFields = fieldnames(segIn.hourly);
    hFields(strcmp(hFields,'time')) = [];

    groups = struct();
    for iF = 1:numel(hFields)
        fname = hFields{iF};
        vec = segIn.hourly.(fname);

        if ~isnumeric(vec)
            continue;
        end

        vec = vec(:);
        if numel(vec) ~= nT
            warning('Skipping %s: length %d != nTime %d', fname, numel(vec), nT);
            continue;
        end

        tokens = regexp(fname,'^(.*)_member(\d+)$','tokens','once');
        if ~isempty(tokens)
            baseVar = tokens{1}; % text before _memberXX
            if ~isfield(groups, baseVar)
                groups.(baseVar) = struct( ...
                    'det', [], ...
                    'memberCols', {{ }}, ...
                    'memberNames', {{ }} );
            end
            groups.(baseVar).memberCols{end+1} = vec;
            groups.(baseVar).memberNames{end+1} = fname;
        else
            % unsuffixed -> deterministic/control/ensemble-mean
            baseVar = fname;
            if ~isfield(groups, baseVar)
                groups.(baseVar) = struct( ...
                    'det', [], ...
                    'memberCols', {{ }}, ...
                    'memberNames', {{ }} );
            end
            groups.(baseVar).det = vec;
        end
    end

    % --- compute stats per base var ---
    stats = struct();
    baseVars = fieldnames(groups);

    for iB = 1:numel(baseVars)
        baseVar = baseVars{iB};
        g = groups.(baseVar);

        if ~isempty(g.memberCols)
            M = cat(2, g.memberCols{:});    % [nT x nMembers]
        else
            M = g.det;                      % [nT x 1]
        end

        if size(M,1) ~= nT
            warning('Skipping var %s due to unexpected matrix shape', baseVar);
            continue;
        end

        nMembers = size(M,2);

        meanVec = mean(M, 2, 'omitnan');
        if nMembers == 1
            stdVec = zeros(nT,1);
            minVec = M;
            maxVec = M;
        else
            stdVec = std(M,0,2,'omitnan');
            minVec = min(M,[],2);
            maxVec = max(M,[],2);
        end

        varStats = struct();
        varStats.raw       = M;          % [nT x nMembers]
        varStats.mean      = meanVec;    % ensemble mean per hour
        varStats.std       = stdVec;     % ensemble spread per hour
        varStats.min       = minVec;     % ensemble min per hour
        varStats.max       = maxVec;     % ensemble max per hour
        varStats.nMembers  = nMembers;   % number of ensemble members
        varStats.det       = g.det;      % deterministic/control line (if present)

        if isstruct(unitsMap) && isfield(unitsMap, baseVar)
            varStats.units = unitsMap.(baseVar);
        else
            varStats.units = '';
        end

        stats.(baseVar) = varStats;
    end

    segOut.stats = stats;
end


%% ========================================================================
function makePlots(seg)
%MAKEPLOTS
% Generate the timeline plots for launch decision making.
% Uses seg.stats.<var>.mean / min / max from processSegment().

    stats = seg.stats;
    tVec  = seg.time;
    tz    = seg.timezone;

    % ---- 1. 10 m wind & gusts ----
    if isfield(stats,'wind_speed_10m')
        figure('Name','Surface Wind (10 m) & Gusts'); clf;
        hold on; grid on; box on;

        ws10 = stats.wind_speed_10m;
        plot(tVec, ws10.mean, 'LineWidth',1.6, 'DisplayName','Wind 10m Mean');
        plot(tVec, ws10.max,  '--', 'LineWidth',1, 'DisplayName','Wind 10m Max');
        plot(tVec, ws10.min,  '--', 'LineWidth',1, 'DisplayName','Wind 10m Min');

        if isfield(stats,'wind_gusts_10m')
            wg10 = stats.wind_gusts_10m;
            plot(tVec, wg10.mean, 'LineWidth',1.6, 'DisplayName','Gust 10m Mean');
        end

        ylabelStr = 'Wind Speed';
        if isfield(ws10,'units') && ~isempty(ws10.units)
            ylabelStr = sprintf('Wind Speed (%s)', ws10.units);
        end
        ylabel(ylabelStr);
        xlabel(sprintf('Time (%s)', tz));
        title('Pad-Level Winds (10 m AGL) and Gusts');
        legend('Location','best');
    end

    % ---- 2. 100 m wind ----
    if isfield(stats,'wind_speed_100m')
        figure('Name','100 m Wind Speed'); clf;
        hold on; grid on; box on;

        ws100 = stats.wind_speed_100m;
        plot(tVec, ws100.mean, 'LineWidth',1.6, 'DisplayName','Wind 100m Mean');
        plot(tVec, ws100.max,  '--', 'LineWidth',1, 'DisplayName','Wind 100m Max');
        plot(tVec, ws100.min,  '--', 'LineWidth',1, 'DisplayName','Wind 100m Min');

        ylabelStr = 'Wind Speed';
        if isfield(ws100,'units') && ~isempty(ws100.units)
            ylabelStr = sprintf('Wind Speed (%s)', ws100.units);
        end
        ylabel(ylabelStr);
        xlabel(sprintf('Time (%s)', tz));
        title('100 m AGL Winds');
        legend('Location','best');
    end

    % ---- 3. Wind direction (10 m) ----
    if isfield(stats,'wind_direction_10m')
        figure('Name','10 m Wind Direction'); clf;
        hold on; grid on; box on;

        wd10 = stats.wind_direction_10m;
        plot(tVec, wd10.mean, 'LineWidth',1.6, 'DisplayName','Dir 10m Mean');

        ylabel('Direction (deg, 0 = N)');
        xlabel(sprintf('Time (%s)', tz));
        title('Wind Direction @10 m AGL');
        legend('Location','best');
    end

    % ---- 4. Pressure / Cloud Cover / Precip ----
    hasPress  = isfield(stats,'surface_pressure');
    hasCloud  = isfield(stats,'cloud_cover');
    hasPrecip = isfield(stats,'precipitation');

    if hasPress || hasCloud || hasPrecip
        figure('Name','Pressure / Clouds / Precip'); clf;
        tl = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

        if hasPress
            nexttile; hold on; grid on; box on;
            sp = stats.surface_pressure;
            plot(tVec, sp.mean, 'LineWidth',1.6);
            ylabelStr = 'Pressure';
            if isfield(sp,'units') && ~isempty(sp.units)
                ylabelStr = sprintf('Pressure (%s)', sp.units);
            end
            ylabel(ylabelStr);
            title('Surface / Sea-Level Pressure');
        end

        if hasCloud
            nexttile; hold on; grid on; box on;
            cc = stats.cloud_cover;
            plot(tVec, cc.mean, 'LineWidth',1.6);
            ylabel('Cloud Cover (%)');
            title('Total Cloud Cover');
        end

        if hasPrecip
            nexttile; hold on; grid on; box on;
            pr = stats.precipitation;
            bar(tVec, pr.mean, 'DisplayName','Precip Mean');
            ylabelStr = 'Precip';
            if isfield(pr,'units') && ~isempty(pr.units)
                ylabelStr = sprintf('Precip (%s)', pr.units);
            end
            ylabel(ylabelStr);
            title('Hourly Precipitation (Ensemble Mean)');
        end

        xlabel(tl, sprintf('Time (%s)', tz));
    end
end