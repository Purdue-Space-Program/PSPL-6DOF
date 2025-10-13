% This script gets the wind data for the input lat and long:
% beware that the namConus data only works for launches in the US.
% Switching to GFS for global launches will be included in the future.

% Wind speed and direction defined by a two-dimensional vector. The
% component u defines the speed of a wind blowing from the West towards the
% East (a negative value therefore implies the opposite direction). The
% component v similarly defines the speed of a wind blowing from the South
% towards the North.

% Units for all other values are given in the struct.

% Input parameters: latitude and longitude
lat = 35.0109889;
lon = -115.4939547;

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

% find the closest time for the 6-DoF for now. In the future this should be
% based on the selection of the user time.

closestTime = weatherData.ts(1);








