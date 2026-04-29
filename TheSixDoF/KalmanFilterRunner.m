data = load("Rocket_IMU_data.mat");

IMU_ES_MEKF(data);

IMU_IEKF(data);

IMU_ES_MEKF_Biases(data);