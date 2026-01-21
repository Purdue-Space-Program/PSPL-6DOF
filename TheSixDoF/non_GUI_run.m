%% %%%%%%%%%% to fix %%%%%%%%%%%%
% rocket.NoseLength % change this to come fron vehicle_parameters.py
% rocket.NoseGeometry % change this to come fron vehicle_parameters.py
% rocket.TotalMass % check if dry or wet?
% launch_site_name % change this to come fron vehicle_parameters.py
% engine.ExitArea = 0.01; change this to come fron vehicle_parameters.py
% engine.ExitPressure = 100000; change this to come fron vehicle_parameters.py
%%%%%%%%%% to fix %%%%%%%%%%%%
%%

vehicle_parameters_csv = fullfile( ...
    "Inputs", ...
    "Saved Rockets", ...
    "FUCK_MATLAB", ...
    "vehicle_parameters.csv");
vehicle_parameters = readtable(vehicle_parameters_csv);

% convert to struct for syntactic sugar
vehicle_parameters = cell2struct( ...
    num2cell(vehicle_parameters.value), ...
    cellstr(vehicle_parameters.parameter_name), ...
    1);


engine = PropulsionSystem("engine_name");

engine.Thrust = vehicle_parameters.jet_thrust;
engine.BurnTime = vehicle_parameters.burn_time;
engine.ExitArea = 0.01;
engine.ExitPressure = 10;

rocket = Rocket("Rocket A");

rocket.TotalLength = vehicle_parameters.total_length; % Vehicle Length [m]
rocket.OuterDiameter = vehicle_parameters.tank_outer_diameter; % Vehicle OD [m]
rocket.NoseLength = 0.5; % Nose Cone Length [m]
rocket.NoseGeometry = "Von Karman"; % Nose Cone Type
rocket.TotalMass = vehicle_parameters.dry_mass; % total mass [kg]
rocket.ComponentList = dictionary(string.empty(0,1), cell.empty(0,1)); % Component Dictionary
% rocket.RASAero_data % RASAero data
rocket.RASAero_data_file_path = "Inputs\Saved Rockets\FUCK_MATLAB\fuck_you6_converted_aero_data.csv"; % The file path of the RASAero data
rocket.CoMOverride = vehicle_parameters.dry_COM_location_from_top; % CoM Override (Dry CoM)
% rocket.CoPOverride % Manual CoP Override


rocket.addComponent(engine)


FAR_latitude = 35.3474; % [degrees]
FAR_longitude = -117.8091; % [degrees]

dairy_farm_latitude = 40.509936707682144; % [degrees]
dairy_farm_longitude = -87.02196060569756; % [degrees]

switch rocket.Name
    case "Rocket A"
        latitude_value = dairy_farm_latitude;
        longitude_value = dairy_farm_longitude;
    case {"CMS", "Copperhead"}
        latitude_value = FAR_latitude;
        longitude_value = FAR_longitude;
    otherwise
        error("Invalid Rocket Name!")
end


launch_site_name = "Launch Trailer";
switch launch_site_name
    case "FAR"
        rail_height = 18.29; % Rail at FAR [m]
    case "Launch Trailer"
        rail_height = 6.096; % Rail on the Launch Trailer [m]
    otherwise
        error("Invalid Launch Site Name!")
end


environment = Environment(latitude_value, longitude_value, datetime("now", "TimeZone","UTC"), rail_height);
environment = getLocalWeather(environment);

settings = IntegratorSettings("apogee", 0.1, "medium");


% selpath = uigetdir(pwd,"Set the root directory to your local 'PSPL-6DOF\TheSixDoF' folder");

Main(rocket, environment, settings, rocket.Name);
% 