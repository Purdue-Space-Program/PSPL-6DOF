% take in the real IMU data:

data = readmatrix("standstill.csv");

IMU_time_ns = data(:,3);

% subtract the first entry to zeroize the time, then convert ns to sec:
IMU_time = (IMU_time_ns - IMU_time_ns(1))/1e9;


IMU_data = struct;

IMU_accel = data(:,4:6);
IMU_gyro = data(:,7:9);
IMU_gps_lla = data(:,10:12);
IMU_gps_v_ned = data(:,13:15);

% find only the switches in GPS
x = data(:,10);  % find indices based on the first column

% Find first non-NaN index
firstValid = find(~isnan(x), 1, 'first');

% Trim leading NaNs
x_valid = x(firstValid:end);

% Find where values change
changeIdx = [1; find(diff(x_valid) ~= 0) + 1];

% Convert back to original indices
changeIdx = changeIdx + firstValid - 1;

% extract all the GPS data at those indices:
gps_lla = IMU_gps_lla(changeIdx, :);
gps_v_ned = IMU_gps_v_ned(changeIdx, :);

gps_time = IMU_time(changeIdx);

% convert the gps lat long alt to ENU

gps_enu = lla2enu(gps_lla,gps_lla(30,:), 'ellipsoid');

% convert velocity to ENU frame:
gps_v_enu(:,1) = gps_v_ned(:,2);
gps_v_enu(:,2) = gps_v_ned(:,1);
gps_v_enu(:,3) = -gps_v_ned(:,3);


% variance and other info:
filename = ['Inputs',filesep,'Saved Rockets',filesep,'CMS.mat'];

rocket = load(filename);
rocket = rocket.rocketObj;

% create a gyroscope, accelerometer, and gps:
gps = GNSS("GPS",1, [4,4,25,0.04,0.04,0.04], .1, 0);

accel = rocket.ComponentList.values{6};
gyro = rocket.ComponentList.values{7};

SamplingRate = 1/100;

gyro.Variance = ones(3,1)*deg2rad(.0035) / sqrt(gyro.SamplingRate); % from VN-200
gyro.Bias0 = 0;
gyro.ScaleFactor = 1.001;
gyro.BiasRandomWalkRate = 2.42e-5;

accel.Variance = 0.0137*ones(3,1);
accel.Bias0 = 0;
accel.ScaleFactor = 1.001;
accel.BiasRandomWalkRate = 3.923e-4;

grav = zeros(numel(time),3);
grav(:,3) = -9.81;



%% Put data into struct

IMU_data.ref_time = IMU_time;
IMU_data.accel_time = IMU_time;
IMU_data.accelometer = IMU_accel;

IMU_data.gyro_time = IMU_time;
IMU_data.gyroscope  = IMU_gyro;

IMU_data.GPS = gps_enu;
IMU_data.GPS_vel = gps_v_enu;
IMU_data.GPS_time = gps_time;

IMU_data.init_cond = [0,0,0,0,0,0,1,0,0,0]'; %pos, vel, quat
IMU_data.grav = grav;

IMU_data.accel_info = accel;
IMU_data.gyro_info = gyro;
IMU_data.gps_info = gps;

IMU_data.ref_traj = zeros(numel(IMU_time),3);

save("Standstill_IMU.mat", "IMU_data")


