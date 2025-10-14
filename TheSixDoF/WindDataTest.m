% This script gets the wind data for the input lat and long:
% beware that the namConus data only works for launches in the US.
% Switching to GFS for global launches will be included in the future.

% Wind speed and direction defined by a two-dimensional vector. The
% component u defines the speed of a wind blowing from the West towards the
% East (a negative value therefore implies the opposite direction). The
% component v similarly defines the speed of a wind blowing from the South
% towards the North.

% Units for all other values are given in the struct.

% Create an environment using the class structure
env = Environment();

% get the location information from the environment
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


%% convert the timestamp to human readable format:

% these dates are milliseconds past Jan 1, 1970
weatherData.ts = weatherData.ts/1000;

weatherData.ts = datetime(weatherData.ts, 'ConvertFrom', 'posixtime');

%% Pull data out of the weatherData list

% find the closest time based on the value set in the environment:

envTime = env.Date;

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

% plot the data
X_values = zeros(1,length(geoHeightInterp));
Y_values = zeros(1,length(geoHeightInterp));
W_values = zeros(1,length(geoHeightInterp));

figure;
subplot(1,2,1)
quiver3(X_values,Y_values,geoHeightInterp, windUInterp,windVInterp,W_values, "off")
xlabel('Wind Vel U [m/s]')
ylabel('Wind Vel V [m/s]')
zlabel('Geopotential Height [m]')
title('Wind Profile')

subplot(1,2,2)
plot(tempInterp,geoHeightInterp)
title('Temperature Profile with Geopotential Height');
xlabel('Temperature [K]')
ylabel('Geopotential Height [m]')
grid on;
