if count(py.sys.path, string(pwd)) == 0
    insert(py.sys.path, int32(0), string(pwd));
end

mod = py.importlib.import_module('OpenMeteoWeather_Request_test');
py.importlib.reload(mod)

weatherData = mod.getData(39.7392,-104.9847)

weatherData = struct(weatherData);