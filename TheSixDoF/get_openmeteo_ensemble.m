function weatherEnsemblePrediction = get_openmeteo_ensemble(lat, lon, dateOrStart, varargin)
%GET_OPENMETEO_ENSEMBLE  Fetch Open-Meteo Ensemble data; store as weatherEnsemblePrediction.
%
% USAGE:
%   % Single day (default vars):
%   weatherEnsemblePrediction = get_openmeteo_ensemble(28.5, -80.6, '2025-11-10');
%
%   % Single day + custom hourly vars (4th arg = vars is OK):
%   weatherEnsemblePrediction = get_openmeteo_ensemble(28.5, -80.6, '2025-11-10', ...
%       {'wind_speed_10m','wind_gusts_10m','wind_speed_100m','wind_direction_100m','pressure_msl','cloud_cover','precipitation','cape'});
%
%   % Single day + placeholder for endDate ([]) + vars:
%   weatherEnsemblePrediction = get_openmeteo_ensemble(28.5, -80.6, '2025-11-10', [], {'wind_speed_10m','wind_speed_100m'});
%
%   % Date RANGE (inclusive) + optional vars:
%   weatherEnsemblePrediction = get_openmeteo_ensemble(28.5, -80.6, '2025-11-10', '2025-11-20', ...
%       {'wind_speed_10m','wind_gusts_10m','wind_speed_100m','wind_direction_100m','pressure_msl','cloud_cover','precipitation','cape'});
%
% DATE WINDOW RULES
%   • Use dates via &start_date= and &end_date=. Same start==end => single day (00:00..23:00 local if timezone=auto).
%   • Max lead ~35 days ahead (model-dependent). Exceeding that results in errors
%
% HOURLY VARIABLES (exact API names):
%   Near-surface: 'wind_speed_10m','wind_gusts_10m','wind_direction_10m','temperature_2m',
%                 'relative_humidity_2m','pressure_msl','cloud_cover','visibility',
%                 'precipitation','rain'
%   Higher Altitude (only 100m tho :( <-- sad face): 'wind_speed_80m','wind_speed_100m','wind_speed_120m',
%                        'wind_direction_80m','wind_direction_100m','wind_direction_120m'
%   Other: 'cape','convective_inhibition','freezing_level_height', 'temperature_850hPa','temperature_500hPa',
%                     'geopotential_height_850hPa','geopotential_height_500hPa'
%
% IMPLEMENTATION
%   - Calls Python script in same folder; fixed model 'gfs_seamless' (long horizon).
%   - Stores result in base workspace: weatherEnsemblePrediction.

    MODELS_FIXED = 'gfs_seamless';

    % Parse optional args (4th can be endDate OR vars)
    endDate = [];
    hourlyVars = [];
    if numel(varargin) == 1
        if iscell(varargin{1}) || (ischar(varargin{1}) && contains(varargin{1}, '_')) || (isstring(varargin{1}) && contains(varargin{1}, "_"))
            hourlyVars = varargin{1};
        else
            endDate = varargin{1};
        end
    elseif numel(varargin) >= 2
        endDate   = varargin{1};
        hourlyVars = varargin{2};
    end

    if isempty(hourlyVars)
        hourlyVars = { ...
            'wind_speed_10m','wind_gusts_10m','wind_direction_10m', ...
            'temperature_2m','relative_humidity_2m','pressure_msl','cloud_cover','visibility','precipitation','rain', ...
            'wind_speed_100m','wind_direction_100m', ...
            'cape','freezing_level_height' ...
        };
    end
    if iscell(hourlyVars)
        hourlyCSV = strjoin(hourlyVars, ',');
    else
        hourlyCSV = char(hourlyVars);
    end

    here = fileparts(mfilename('fullpath'));
    pyscript = fullfile(here, 'openmeteo_ensemble.py');
    if ~isfile(pyscript)
        error('Python script not found at: %s', pyscript);
    end

    if isempty(endDate)
        cmd = sprintf('python "%s" --lat %g --lon %g --date "%s" --models "%s" --hourly "%s" --timezone auto', ...
            pyscript, lat, lon, dateOrStart, MODELS_FIXED, hourlyCSV);
    else
        cmd = sprintf('python "%s" --lat %g --lon %g --start-date "%s" --end-date "%s" --models "%s" --hourly "%s" --timezone auto', ...
            pyscript, lat, lon, dateOrStart, endDate, MODELS_FIXED, hourlyCSV);
    end

    [status, out] = system(cmd);
    if status ~= 0
        error('Open-Meteo call failed:\n%s', out);
    end

    weatherEnsemblePrediction = jsondecode(out); 
    assignin('base', 'weatherEnsemblePrediction', weatherEnsemblePrediction);
end
