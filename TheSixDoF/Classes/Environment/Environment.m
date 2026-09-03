classdef Environment
    % The environment class models the environment system for the rocket,
    % including the sea level height of the location, local gravity, and
    % atmospheric conditions at the current time.

    properties
        Date (1,1) datetime = datetime("now", "TimeZone", "UTC");
        Elevation (1,1) double = 627.91;
        LatLong (1,2) double = [35.347444074690735, -117.8090720168799]
        railHeight (1,1) double = 18.29; % FAR rail height [m]
        Atmosphere (:,5) double
        Wind (:,3) double
        LaunchWeather struct = struct();   % ADDED: full Open-Meteo forecast w/ ensemble stats & times
    end

    methods
        function env = Environment(lat, long, date)
            % Environment creates a new launch site given an input
            % environment and the position in latitude and longitude. The
            % function automatically generates the elevation of the
            % location.
            arguments
                lat (1,1) double = 35.347444074690735
                long (1,1) double = -117.8090720168799
                date (1,1) datetime = datetime("now", "TimeZone", "UTC")
            end
            env.LatLong = [lat,long];
            try
                env.Elevation = getElevation(env);
            catch ME
                warning('Environment:elevationFailed', ...
                    'getElevation failed (%s). Set env.Elevation manually.', ME.message);
                env.Elevation = 0;
            end
            env.Date = date;
        end

        function elev = getElevation(env)
            % the open-meteo also gives elevation, may switch over to this.
            loc = txsite("Latitude", env.LatLong(1),"Longitude",env.LatLong(2));
            elev = elevation(loc);
        end

        function env = getLocalWeather(env)
            % get the local weather based on the environment that has been
            % set by the user.
            %
            % When called with 1 outputs, n x 3 vector containing the
            % geopotential height, wind in the u direction, and wind in the
            % v direction.
            %
            % When called with 2 outputs, the first variable returns n x 1
            % vector of the geopotential height, and the 2nd returns n x 2
            % vector of the wind in the u and v direction.
            %
            % With three to five outputs, the output looks like:
            % [geoHeight, wind, temp, surfGust, surfPressure]

            % Wind speed and direction defined by a two-dimensional vector. The
            % component u defines the speed of a wind blowing from the West towards the
            % East (a negative value therefore implies the opposite direction). The
            % component v similarly defines the speed of a wind blowing from the South
            % towards the North.

            % Units for all other values are given in the struct.

            % get the location information from the environment
            arguments
                env Environment
            end

            latlong = env.LatLong;

            lat = latlong(1);
            lon = latlong(2);

            % get the folder. Would be better to do this dynamically
            folder = fullfile(pwd, 'TheSixDoF');

            if count(py.sys.path, string(pwd)) == 0
                insert(py.sys.path, int32(0), string(folder));
            end

            % import the python file which contains all of the scripts:
            mod = py.importlib.import_module('OpenMeteoWeatherRequest');
            py.importlib.reload(mod)

            hourlyWeatherData = mod.getHourlyWeather(lat,lon);

            hourlyWeatherData = struct(hourlyWeatherData);

            % Convert all Python types to native MATLAB types
            f = fieldnames(hourlyWeatherData);
            for i = 1:numel(f)
                val = hourlyWeatherData.(f{i});
                if isa(val, 'py.numpy.ndarray')
                    hourlyWeatherData.(f{i}) = double(val);
                elseif isa(val, 'py.list')
                    try
                        hourlyWeatherData.(f{i}) = double(py.array.array('d', val));
                    catch
                        % list has non-numeric elements — convert element-wise
                        c = cell(val);
                        hourlyWeatherData.(f{i}) = cellfun(@(x) double(x), c);
                    end
                elseif isa(val, 'py.float') || isa(val, 'py.int')
                    hourlyWeatherData.(f{i}) = double(val);
                end
                % py.str and other non-numeric types: leave as-is
            end

            hourlyWeatherData.date = datetime(hourlyWeatherData.date, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');

            % convert to local time:
            hourlyWeatherData.dateLocal = hourlyWeatherData.date + seconds(hourlyWeatherData.timeOffset);

            % update the local environment time
            envTime = env.Date;

            % Find the index of the closest timestamp in weatherData
            [~, closestIndex] = min(abs(hourlyWeatherData.date - envTime));

            closestTime = hourlyWeatherData.date(closestIndex);

            % get the single variable fields:
            surfPres = hourlyWeatherData.surface_pressure(closestIndex);

            surfGust = hourlyWeatherData.wind_gusts_10m(closestIndex);

            % Extract each pressure level explicitly by constructing field names.
            % This avoids two bugs in the field-search approach:
            %   1) MATLAB struct fields are alphabetically ordered, so
            %      contains()-based loops don't match presLevels order.
            %   2) contains(fields,'temperature_') also matches temperature_2m,
            %      giving one extra element and breaking all downstream indexing.
            presLevelsHPa = [1000, 925, 850, 700, 500, 300, 200, 100, 50];
            nLev = numel(presLevelsHPa);

            geoHeight_raw  = zeros(1, nLev);
            windSpeed_raw  = zeros(1, nLev);
            windDir_raw    = zeros(1, nLev);
            tempKelvin_raw = zeros(1, nLev);

            for idx = 1:nLev
                p = presLevelsHPa(idx);
                geoHeight_raw(idx)  = hourlyWeatherData.(sprintf('geopotential_height_%dhPa', p))(closestIndex);
                windSpeed_raw(idx)  = hourlyWeatherData.(sprintf('wind_speed_%dhPa',          p))(closestIndex);
                windDir_raw(idx)    = hourlyWeatherData.(sprintf('wind_direction_%dhPa',       p))(closestIndex);
                tempKelvin_raw(idx) = hourlyWeatherData.(sprintf('temperature_%dhPa',          p))(closestIndex) + 273.15;
            end

            % Keep only levels at or above site elevation, sorted ascending
            keep = geoHeight_raw >= env.Elevation;
            geoHeight  = geoHeight_raw(keep);
            windSpeed  = windSpeed_raw(keep);
            windDir    = windDir_raw(keep);
            tempKelvin = tempKelvin_raw(keep);
            presLevels = (presLevelsHPa(keep) * 100)';

            [geoHeight, sIdx] = sort(geoHeight);
            windSpeed  = windSpeed(sIdx);
            windDir    = windDir(sIdx);
            tempKelvin = tempKelvin(sIdx);
            presLevels = presLevels(sIdx);

            % interpolate the data for more fine grain height increments:
            geoHeightInterp = linspace(geoHeight(1),geoHeight(end),1000);

            % interpolate the other data:
            % Interpolate wind data u and v using the same height interpolation
            windSpeedInterp = interp1(geoHeight, windSpeed, geoHeightInterp, 'pchip', 'extrap');
            windDirInterp = interp1(geoHeight, windDir, geoHeightInterp, 'pchip', 'extrap');
            tempInterp = interp1(geoHeight, tempKelvin, geoHeightInterp, 'linear', 'extrap');
            presInterp = interp1(geoHeight, presLevels, geoHeightInterp, 'pchip', 'extrap');
            rhoInterp = presInterp ./ (constant.R_AIR * tempInterp);
            aInterp = sqrt(constant.R_AIR .* constant.GAMMA_AIR .* tempInterp);

            % get the wind directions. Get it as the direction the wind is actually
            % going (not the direction it is coming from)!
            % U is the N-S direction. A positive value is wind blowing from
            % the south (the wind is going north)
            % V is the E-W direction. A positive value is wing blowing from
            % the west (going east)

            windU = -windSpeedInterp .* cosd(windDirInterp);
            windV = -windSpeedInterp .* sind(windDirInterp);

            % put this into the environment:
            env.Atmosphere = [geoHeightInterp', tempInterp', aInterp', presInterp', rhoInterp'];

            env.Wind = [geoHeightInterp', windU', windV'];

            % change the number of outputs based on the users choice.
            % Helpful if the user wants outputs for something else:
            % nOutputs = nargout;

            % switch nOutputs
            % 
            %     case 1
            %         varargout{1} = [geoHeightInterp;windU;windV]';
            %     case 2
            %         varargout{1} = geoHeightInterp';
            %         varargout{2} = [windU;windV]';
            %     case 3
            %         varargout{1} = geoHeightInterp';
            %         varargout{2} = [windU;windV]';
            %         varargout{3} = tempInterp';
            %     case 4
            %         varargout{1} = geoHeightInterp';
            %         varargout{2} = [windUInterp;windVInterp]';
            %         varargout{3} = tempInterp';
            %         varargout{4} = surfGust;
            %     case 5
            %         varargout{1} = geoHeightInterp';
            %         varargout{2} = [windUInterp;windVInterp]';
            %         varargout{3} = tempInterp';
            %         varargout{4} = surfGust;
            %         varargout{5} = surfPres;
            % end

        end

        function env = getLaunchWeather(env, startDate, endDate, hourlyVars)
            % getLaunchWeather
            % (MORE INFO IN get_openmeteo_ensemble.m)
            % Pull Open-Meteo launch weather/forecast (or historical, if in the past)
            % and store it on env.LaunchWeather.
            %
            % After calling this method:
            %   env.LaunchWeather
            %       .ensemble (for today/future) or .archive (for past)
            %       .ensemble.time  -> datetime array for each forecast hour
            %       .ensemble.stats.wind_speed_10m.mean
            %       .ensemble.stats.wind_gusts_10m.max
            %       .ensemble.stats.wind_speed_100m.mean
            %       .ensemble.stats.cloud_cover.mean
            %       .ensemble.stats.precipitation.mean
            %
            % NOTE: get_openmeteo_ensemble() also drops a copy into
            % base workspace as weatherEnsemblePrediction for convenience.
            %
            % EXAMPLE:
            %   env = env.getLaunchWeather();  % use env.Date
            %   env = env.getLaunchWeather('2025-11-10');         % that single day
            %   env = env.getLaunchWeather('2025-11-10','2025-11-12'); % 3-day window

            % Default argument handling
            if nargin < 2 || isempty(startDate)
            % default: use env.Date rounded to date string
            startDate = datestr(env.Date, 'yyyy-mm-dd');
            end
    
            if nargin < 3 || isempty(endDate)
                endDate = startDate;
            end
    
            if nargin < 4
                hourlyVars = [];
            end
    
            % Call the main fetcher/plotter
            wx = get_openmeteo_ensemble( ...
                    env.LatLong(1), ...   % lat
                    env.LatLong(2), ...   % lon
                    startDate, ...
                    endDate, ...
                    hourlyVars);
    
            % Stash the full struct
            env.LaunchWeather = wx;
    
            % Also mirror the pieces we care about for sim:
            if isfield(wx,'Atmosphere')
                env.Atmosphere = wx.Atmosphere;
            end
            if isfield(wx,'Wind')
                env.Wind = wx.Wind;
            end
        end
    end
end