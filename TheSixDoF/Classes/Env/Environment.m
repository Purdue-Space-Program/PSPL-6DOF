classdef Environment
    % The environment class models the environment system for the rocket,
    % including the sea level height of the location, local gravity, and

    properties
        Date (1,1) datetime = datetime("now");
        Elevation (1,1) double = 627.91;
        LatLong (1,2) double = [35.347444074690735, -117.8090720168799]
        geocentricRadius (1,1) double = 6.371077849286893e6;
        railHeight (1,1) double = 18.29; % FAR rail height [m]
    end

    methods
        function env = Environment(lat, long, date)
            % setLaunchSite creates a new launch site given an input
            % environment and the position in latitude and longitude. The
            % function automatically generates the elevation
                env.LatLong = [lat,long];
                env.Elevation = getElevation(env);
                env.Date = date;
        end

        function elev = getElevation(env)
            loc = txsite("Latitude", env.LatLong(1),"Longitude",env.LatLong(2));
            elev = elevation(loc);
        end

        function varargout = getLocalWeather(env)
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

            % Path to your Python script
            python_script = 'WindyAPI_Request.py';  % Adjust to the actual location

            % Construct system command
            cmd = sprintf('python "%s" %.6f %.6f', python_script, lat, lon);

            % Run the Python script and capture output
            [status, result] = system(cmd);

            % Check result
            if status == 0
                % Clean up whitespace/newlines and decode JSON
                result = strtrim(result);
                weatherData = jsondecode(result);
            else
                fprintf('Error running Python script:\n%s\n', result);
            end


            % convert the timestamp to human readable format and make it a datetime object:
            % these dates are milliseconds past Jan 1, 1970. Convert these
            % with the datetime function.
            weatherData.ts = weatherData.ts/1000;

            weatherData.ts = datetime(weatherData.ts, 'ConvertFrom', 'posixtime');

            % Pull data out of the weatherData list
            % find the closest time based on the value set in the environment:
            envTime = env.Date;

            % if the requested date is more than 1 day in the past
            if env.Date - datetime("today") < 0.5
                warndlg('Pulling historic weather data is currently not supported. The script will pull the most current weather data.')
                warning('Pulling historic weather data is currently not supported. The script will pull the most current weather data.')
            end

            % Find the index of the closest timestamp in weatherData
            [~, closestIndex] = min(abs(weatherData.ts - envTime));

            closestTime = weatherData.ts(closestIndex);

            % Get all field names in your struct
            fields = fieldnames(weatherData);

            % Keep only the geopotential height fields
            ghFields = fields(contains(fields, 'gh_'));

            geoHeight(1) = env.Elevation;

            % Loop and extract the value at that index
            for idx = 2:numel(ghFields)
                f = ghFields{idx};
                geoHeight(idx) = weatherData.(f)(closestIndex);  % access by dynamic field name
            end

            FieldsFilter = geoHeight >= env.Elevation;


            % remove any of the geoheights which are less than the site elevation
            geoHeight(geoHeight < env.Elevation) = [];

            % use the fields filter to go through the rest of the data and pull the
            % appropriate values:

            % wind data u (west towards east):
            windUFields = fields(contains(fields, 'wind_u'));

            windUFields = windUFields(FieldsFilter);

            % Loop and extract the value at that index
            for idx = 1:numel(windUFields)
                f = windUFields{idx};
                windU(idx) = weatherData.(f)(closestIndex);  % access by dynamic field name
            end

            % wind data v (south towards north):
            windVFields = fields(contains(fields, 'wind_v'));

            windVFields = windVFields(FieldsFilter);

            % Loop and extract the value at that index
            for idx = 1:numel(windVFields)
                f = windVFields{idx};
                windV(idx) = weatherData.(f)(closestIndex);  % access by dynamic field name
            end

            % temp data

            tempFields = fields(contains(fields, 'temp_'));

            tempFields = tempFields(FieldsFilter);

            % Loop and extract the value at that index
            for idx = 1:numel(tempFields)
                f = tempFields{idx};
                tempKelvin(idx) = weatherData.(f)(closestIndex);  % access by dynamic field name
            end

            % interpolate the data for more data
            geoHeightInterp = linspace(geoHeight(1),geoHeight(end),100);

            % interpolate the other data:
            % Interpolate wind data u and v using the same height interpolation
            windUInterp = interp1(geoHeight, windU, geoHeightInterp, 'linear', 'extrap');
            windVInterp = interp1(geoHeight, windV, geoHeightInterp, 'linear', 'extrap');
            tempInterp = interp1(geoHeight, tempKelvin, geoHeightInterp, 'linear', 'extrap');

            % change the number of outputs based on the user input:
            nOutputs = nargout;

            switch nOutputs

                case 1
                    varargout{1} = [geoHeightInterp;windUInterp;windVInterp]';
                case 2
                    varargout{1} = geoHeightInterp';
                    varargout{2} = [windUInterp;windVInterp]'; 
                case 3
                    varargout{1} = geoHeightInterp';
                    varargout{2} = [windUInterp;windVInterp]'; 
                    varargout{3} = tempInterp';
                case 4
                    gust = weatherData.gust_surface(closestIndex);

                    varargout{1} = geoHeightInterp';
                    varargout{2} = [windUInterp;windVInterp]'; 
                    varargout{3} = tempInterp';
                    varargout{4} = gust;
                case 5
                    gust = weatherData.gust_surface(closestIndex);
                    pres = weatherData.pressure_surface(closestIndex);

                    varargout{1} = geoHeightInterp';
                    varargout{2} = [windUInterp;windVInterp]'; 
                    varargout{3} = tempInterp';
                    varargout{4} = gust;
                    varargout{5} = pres;
            end

        end
    end
end