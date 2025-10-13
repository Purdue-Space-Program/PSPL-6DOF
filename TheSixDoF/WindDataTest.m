% Input parameters: latitude and longitude
lat = input('Enter latitude: ');
lon = input('Enter longitude: ');

% Set the path to your Python script (replace with the actual path)
python_script = 'WindyAPI_Request.py';  % For example, 'C:/my_scripts/fetch_weather.py'

% Run the Python script with pyrunfile, and capture the printed output
[status, result] = system(['python ', python_script, ' ', num2str(lat), ' ', num2str(lon)]);

% If the Python script executed successfully, result should contain the JSON output
if status == 0
    % Parse the result (JSON string) into MATLAB data
    weather_data = jsondecode(result);
    
    % Display the result (weather data)
    disp(weather_data);
else
    disp('Error in running the Python script:');
    disp(result);
end
