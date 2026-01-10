if count(py.sys.path, string(pwd)) == 0
    insert(py.sys.path, int32(0), string(pwd));
end

% import the python file which contains all of the scripts:
mod = py.importlib.import_module('OpenMeteoWeatherRequest');
py.importlib.reload(mod);

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

% update the local environment time
envTime = env.Date;

% Find the index of the closest timestamp in weatherData
[~, closestIndex] = min(abs(hourlyWeatherData.date - envTime));

closestTime = hourlyWeatherData.date(closestIndex);

% Get all field names in the struct
fields = fieldnames(hourlyWeatherData);

% Keep only the geopotential height fields
ghFields = fields(contains(fields, 'geopotential_height_'));

geoHeight(1) = env.Elevation;

% Loop and extract the value at that index
for idx = 1:numel(ghFields)
    f = ghFields{idx};
    geoHeight(idx) = hourlyWeatherData.(f)(closestIndex);
end

PresLevels = [1000e2, 975e2, 950e2, 925e2, 900e2, 850e2, 800e2, 700e2, 600e2, ...
    500e2, 400e2, 300e2, 250e2, 200e2, 150e2, 100e2, 70e2, 50e2, 30e2]';

FieldsFilter = geoHeight >= env.Elevation;

PresLevels = PresLevels(FieldsFilter);

% remove any of the geoheights which are less than the site elevation
geoHeight(geoHeight < env.Elevation) = [];

% use the fields filter to go through the rest of the data and pull the
% appropriate values:

% wind speed data:
windSpeedFields = fields(contains(fields, 'wind_speed_'));

windSpeedFields = windSpeedFields(FieldsFilter);

% Loop and extract the value at that index
for idx = 1:numel(windSpeedFields)
    f = windSpeedFields{idx};
    windSpeed(idx) = hourlyWeatherData.(f)(closestIndex);
end

% wind dir data:
windDirFields = fields(contains(fields, 'wind_direction_'));

windDirFields = windDirFields(FieldsFilter);

% Loop and extract the value at that index
for idx = 1:numel(windDirFields)
    f = windDirFields{idx};
    windDir(idx) = hourlyWeatherData.(f)(closestIndex);
end

% temp data

tempFields = fields(contains(fields, 'temperature_'));

tempFields = tempFields(FieldsFilter);

% Loop and extract the value at that index
for idx = 1:numel(tempFields)
    f = tempFields{idx};
    tempKelvin(idx) = hourlyWeatherData.(f)(closestIndex) + 273.15; 
end

% interpolate the data for more data
geoHeightInterp = linspace(geoHeight(1),geoHeight(end),300);

% interpolate the other data:
% Interpolate wind data u and v using the same height interpolation
windSpeedInterp = interp1(geoHeight, windSpeed, geoHeightInterp, 'pchip', 'extrap');
windDirInterp = interp1(geoHeight, windDir, geoHeightInterp, 'pchip', 'extrap');
tempInterp = interp1(geoHeight, tempKelvin, geoHeightInterp, 'pchip', 'extrap');
presInterp = interp1(geoHeight, PresLevels, geoHeightInterp, 'pchip', 'extrap');
R = 287.05;
rhoInterp = presInterp ./ (R * tempInterp);

% plot the data
X_values = zeros(1,length(geoHeightInterp));
Y_values = zeros(1,length(geoHeightInterp));
W_values = zeros(1,length(geoHeightInterp));

% get the wind directions. Get it as the direction the wind is actually
% going (not the direction it is coming from)!
windU = windSpeedInterp .* cosd(windDirInterp);
windV = windSpeedInterp .* sind(windDirInterp);

figure;
subplot(1,2,1)
quiver3(X_values,Y_values,geoHeightInterp,windU,windV,W_values, "off")
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
hold on

% compare the results against standard atmosphere:
[T, a, P, rho] = atmosisa(geoHeightInterp);

figure;
subplot(1,3,1)
plot(tempInterp,geoHeightInterp)
hold on
plot(T, geoHeightInterp)
title('Temperature Profile');
xlabel('Temperature [K]')
ylabel('Geopotential Height [m]')
grid on;
legend('API Data', 'Standard Atmosphere')


% Compare the pressure levels against standard atmosphere
subplot(1,3,2)
plot(presInterp, geoHeightInterp)
hold on
plot(P, geoHeightInterp)
title('Pressure Profile');
xlabel('Pressure [Pa]')
ylabel('Geopotential Height [m]')
grid on;
legend('API Data', 'Standard Atmosphere')

% compare the density levels against standard atmosphere
subplot(1,3,3)
plot(rhoInterp, geoHeightInterp)
hold on
plot(rho, geoHeightInterp)
title('Density Profile');
xlabel('Density [$kg/m^3$]')
ylabel('Geopotential Height [m]')
grid on;
legend('API Data', 'Standard Atmosphere')







