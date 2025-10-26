if count(py.sys.path, string(pwd)) == 0
    insert(py.sys.path, int32(0), string(pwd));
end

% import the python file which contains all of the scripts:
mod = py.importlib.import_module('OpenMeteoWeatherRequest');
py.importlib.reload(mod)

latlong = env.LatLong;

lat = latlong(1);
lon = latlong(2);

hourlyWeatherData = mod.getHourlyWeather(lat,lon);

currentWeatherData = mod.getCurrentWeather(lat,lon);

hourlyWeatherData = struct(hourlyWeatherData);
currentWeatherData = struct(currentWeatherData);

% adjust the current weather:
currentWeatherData.time = datetime("now", "TimeZone", "UTC");

% Auto-convert all py.numpy arrays or lists to MATLAB doubles
f = fieldnames(hourlyWeatherData);
for i = 1:numel(f)
    val = hourlyWeatherData.(f{i});
    if isa(val, 'py.numpy.ndarray')
        hourlyWeatherData.(f{i}) = double(val);
    elseif isa(val, 'py.list')
        hourlyWeatherData.(f{i}) = double(py.array.array('d', val));
    end
end

hourlyWeatherData.date = datetime(hourlyWeatherData.date, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');

% convert to local time:
hourlyWeatherData.dateLocal = hourlyWeatherData.date + seconds(hourlyWeatherData.timeOffset);

% get the data from the environment:
envTime = env.Date;

% Find the index of the closest timestamp in weatherData
[~, closestIndex] = min(abs(hourlyWeatherData.date - envTime));

closestTime = hourlyWeatherData.date(closestIndex);

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




