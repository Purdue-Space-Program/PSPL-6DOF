 classdef constant
    properties (Constant = true)

         % conversion factors
         mu = 3.986004418e14;    % gravitional constant of earth [m^3/s^-2]
         m_ft = 3.28084;         % meters to feet conversion
         ft_m = 0.3048;          % ft to m conversion
         m_inch = 39.37;         % m to in conversion
         lbm_kg = 0.45359237;    % lbm to kg conversion
         in3_m3 = 0.00001639;    % in3 to m3 conversion

         % properties of gases:
         R_u = 8.31446261815324; % universal gas constant [J·K/mol]


         % Gas properties for air:
         MW_air = 28.9647/1000;  % weight of air [kg/mol]

         % specific gas constant for dry air [J/(kg·K)]
         R_air = constant.R_u / constant.MW_air;
         gamma_air = 1.4;

    end
 end