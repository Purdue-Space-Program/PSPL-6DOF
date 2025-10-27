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
% OUTPUT STRUCT: weatherEnsemblePrediction
%   .raw       : the exact JSON from Python (archive / ensemble / or both)
%   .archive   : processed past segment (if any)
%   .ensemble  : processed future segment (if any)
%   .plot_info : which segment we plotted + which vars we asked for
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
% SIDE EFFECTS:
%   - Generates figures for winds/gusts/shear/clouds/rain/pressure
%   - Assigns the final struct into base workspace as weatherEnsemblePrediction
%
% REQUIREMENTS:
%   - Python 3 on PATH with "requests" installed
%   - openmeteo_ensemble.py in the same folder as this .m file

    % --------------------------
    % Arg handling / defaults
    % --------------------------
    if nargin < 3
        error('Usage: get_openmeteo_ensemble(lat, lon, startDate [,endDate] [,hourlyVars])');
    end

    if nargin < 4 || isempty(endDate)
        endDate = startDate;
    end

    defaultHourlyVars = { ...
        'wind_speed_10m', ...      % m/s at 10 m AGL (pad winds, launch commit)
        'wind_speed_100m', ...     % m/s at 100 m AGL (low-level shear / max-Q)
        'wind_gusts_10m', ...      % gust @10 m (peak last hr; critical for pad ops)
        'wind_direction_10m', ...  % deg true at 10 m (crosswind / plume drift)
        'surface_pressure', ...    % hPa near surface / MSL pressure (fronts)
        'cloud_cover', ...         % % sky obscured (ceiling / lightning rules)
        'precipitation' ...        % mm/hr total precip (no rain near pad)
    };

    if nargin < 5 || isempty(hourlyVars)
        hourlyVars = defaultHourlyVars;
    end
    if ischar(hourlyVars)
        hourlyVars = {hourlyVars};
    end

    % We'll pass hourly vars to Python as separate tokens.
    hourlyArgStr = strjoin(hourlyVars, ' ');


    % Locate python helper
    thisFilePath = mfilename('fullpath');
    [thisDir,~,~] = fileparts(thisFilePath);
    pyScriptPath = fullfile(thisDir,'openmeteo_ensemble.py');

    if ~isfile(pyScriptPath)
        error('Python helper not found at: %s', pyScriptPath);
    end

    % Ensemble model for FUTURE side. We pick a long-range GFS-style ensemble,
    % which can extend up to ~35 days ahead at lower resolution. :contentReference[oaicite:6]{index=6}
    modelName = 'gfs_seamless';

    % Build python command. We always pass both start and end.
    cmd = sprintf([ ...
        'python "%s" --lat %g --lon %g ' , ...
        '--start-date "%s" --end-date "%s" ' , ...
        '--models "%s" ' , ...
        '--timezone auto ' , ...
        '--wind-speed-unit ms ' , ...         % m/s is nice for structural limits
        '--precipitation-unit mm ' , ...      % mm/hr
        '--temperature-unit celsius ' , ...
        '--hourly %s' ], ...
        pyScriptPath, lat, lon, startDate, endDate, modelName, hourlyArgStr);

    % Call Python, decode JSON
    [status, cmdout] = system(cmd);
    if status ~= 0
        error('Open-Meteo call failed:\n%s', cmdout);
    end

    try
        rawData = jsondecode(cmdout);
    catch decodeErr
        error('Failed to jsondecode python output.\nRaw output:\n%s\n\nError:\n%s', ...
              cmdout, decodeErr.message);
    end

    weatherEnsemblePrediction = struct();
    weatherEnsemblePrediction.raw = rawData;

    % The Python can return:
    %   {archive: <seg>, ensemble: <seg>}  OR just <seg>.
    hasArchive = isstruct(rawData) && isfield(rawData,'archive');
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

    if isempty(weatherEnsemblePrediction.archive) && isempty(weatherEnsemblePrediction.ensemble)
        % If python already "unwrapped" because there's only one segment
        if isstruct(rawData) && isfield(rawData,'hourly')
            weatherEnsemblePrediction.ensemble = processSegment(rawData);
        end
    end

    % Pick which chunk to plot:
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

    % Plots (winds, gusts, shear, direction, pressure/clouds/rain)
    if ~isempty(plotSeg)
        makePlots(plotSeg);
    else
        warning('No data segment available to plot.');
    end

    % Make it easy to poke around interactively:
    assignin('base','weatherEnsemblePrediction', weatherEnsemblePrediction);
end


function segOut = processSegment(segIn)

    segOut = struct();

    if ~isfield(segIn,'hourly')
        segOut = [];
        return;
    end

    % --- Pull & convert time ---
    timeStrings = segIn.hourly.time;
    if iscell(timeStrings), timeStrings = string(timeStrings);
    elseif ischar(timeStrings), timeStrings = string(timeStrings);
    else, timeStrings = string(timeStrings);
    end

    % timezone from API (e.g. 'America/New_York')
    tz = 'UTC';
    if isfield(segIn,'timezone')
        if ischar(segIn.timezone)
            tz = segIn.timezone;
        elseif isstring(segIn.timezone)
            tz = char(segIn.timezone);
        end
    end

    % Parse timestamps as local time in that tz.
    % Open-Meteo returns "YYYY-MM-DDTHH:MM"
    try
        tVec = datetime(timeStrings, ...
                        'InputFormat',"yyyy-MM-dd'T'HH:mm", ...
                        'TimeZone',tz);
    catch
        tVec = datetime(timeStrings,'TimeZone',tz);
    end
    tVec = tVec(:); % column
    nT = numel(tVec);

    segOut.time = tVec;
    segOut.timezone = tz;

    % --- Units map if provided ---
    unitsMap = struct();
    if isfield(segIn,'hourly_units')
        unitsMap = segIn.hourly_units;
    end

    % --- Group members by base variable name ---
    hFields = fieldnames(segIn.hourly);
    hFields(strcmp(hFields,'time')) = [];

    % We'll build a struct of groups:
    % groups.(baseVar).det         -> [nT x 1] (control/mean series w/out _memberXX)
    % groups.(baseVar).memberCols  -> { [nT x 1], [nT x 1], ... }
    % groups.(baseVar).memberNames -> {'wind_speed_10m_member01', ...}
    groups = struct();

    for iF = 1:numel(hFields)
        fname = hFields{iF};
        vec = segIn.hourly.(fname);

        % force column vector double
        if ~isnumeric(vec)
            % skip weird/non-numeric series
            continue;
        end
        vec = vec(:);

        if numel(vec) ~= nT
            % dimension mismatch is suspicious: skip
            warning('Skipping %s: length %d does not match nTime=%d', ...
                fname, numel(vec), nT);
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
            % unsuffixed -> treat as deterministic/control/ensemble-mean line
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

    % For each baseVar, build matrix and compute per-hour stats across members
    stats = struct();
    baseVars = fieldnames(groups);

    for iB = 1:numel(baseVars)
        baseVar = baseVars{iB};
        g = groups.(baseVar);

        % Stack all member columns into [nT x nMembers]
        if ~isempty(g.memberCols)
            M = cat(2, g.memberCols{:});    % each col = one ensemble member track
        else
            % no members? fall back to deterministic line alone
            M = g.det;
        end

        % Safety: ensure row count matches nT
        if size(M,1) ~= nT
            warning('Skipping var %s due to unexpected matrix shape', baseVar);
            continue;
        end

        nMembers = size(M,2);

        % Per-hour stats across ensemble members (columns)
        % This is what we want for launch decision support:
        % at each hour, what does the ensemble say? mean / spread / extremes.
        meanVec = mean(M, 2, 'omitnan');     % [nT x 1]
        if nMembers == 1
            stdVec = zeros(nT,1);
            minVec = M;
            maxVec = M;
        else
            stdVec = std(M, 0, 2, 'omitnan');% [nT x 1]
            minVec = min(M, [], 2);          % [nT x 1]
            maxVec = max(M, [], 2);          % [nT x 1]
        end

        varStats = struct();
        varStats.raw       = M;              % [nT x nMembers], all ensemble members
        varStats.mean      = meanVec;        % ensemble mean per hour
        varStats.std       = stdVec;         % ensemble spread per hour
        varStats.min       = minVec;         % ensemble min per hour
        varStats.max       = maxVec;         % ensemble max per hour
        varStats.nMembers  = nMembers;       % how many ensemble members we actually got (should be 31)
        varStats.det       = g.det;          % "control"/unsuffixed series for reference

        % Attach units if we have them (unitsMap usually keys on baseVar like wind_speed_10m)
        if isstruct(unitsMap) && isfield(unitsMap, baseVar)
            varStats.units = unitsMap.(baseVar);
        else
            varStats.units = '';
        end

        stats.(baseVar) = varStats;
    end

    segOut.stats = stats;
end


function makePlots(seg)
%MAKEPLOTS Quick-look launch weather plots:
%   1. 10 m wind & gusts
%   2. 100 m wind (low-level shear)
%   3. 10 m wind direction
%   4. Pressure / cloud cover / precip
%
% We use the per-hour ensemble stats we just computed:
%   .mean = mean across members
%   .min/.max = envelope across members
%   .std = spread (uncertainty)

    stats = seg.stats;
    tVec  = seg.time;
    tz    = seg.timezone;

    % -------- 1. Surface wind & gusts (10 m AGL) --------
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

    % -------- 2. 100 m wind speed (ascent shear environment) --------
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
        title('100 m AGL Winds (Low-Level Shear / Max-Q)');
        legend('Location','best');
    end

    % -------- 3. 10 m wind direction --------
    if isfield(stats,'wind_direction_10m')
        figure('Name','10 m Wind Direction'); clf;
        hold on; grid on; box on;

        wd10 = stats.wind_direction_10m;
        plot(tVec, wd10.mean, 'LineWidth',1.6, 'DisplayName','Dir 10m Mean');
        % direction spread (std) is less meaningful past wrap-around,
        % but you can still eyeball mean trend for plume drift.
        ylabel('Direction (deg, 0 = N)');
        xlabel(sprintf('Time (%s)', tz));
        title('Wind Direction @10 m AGL');
        legend('Location','best');
    end

    % -------- 4. Pressure / Clouds / Precip --------
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
