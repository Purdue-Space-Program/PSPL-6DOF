if count(py.sys.path, string(pwd)) == 0
    insert(py.sys.path, int32(0), string(pwd));
end

% import the python file which contains all of the scripts:
mod = py.importlib.import_module('OpenMeteoWeatherRequest');
py.importlib.reload(mod)

% if the user requests the current launch data, run the 

lat = 39.7392;
lon = -104.9847;

weatherData = mod.getWeather(lat,lon);

weatherData = struct(weatherData);

% Auto-convert all py.numpy arrays or lists to MATLAB doubles
f = fieldnames(weatherData);
for i = 1:numel(f)
    val = weatherData.(f{i});
    if isa(val, 'py.numpy.ndarray')
        weatherData.(f{i}) = double(val);
    elseif isa(val, 'py.list')
        weatherData.(f{i}) = double(py.array.array('d', val));
    end
end

weatherData.date = datetime(weatherData.date, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');

% convert to local time:
weatherData.dateLocal = weatherData.date + seconds(weatherData.timeOffset);
