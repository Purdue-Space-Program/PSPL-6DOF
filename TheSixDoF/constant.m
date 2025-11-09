 classdef constant
    properties (Constant = true)

         % Conversion Factors
         M_TO_FT = 3.28084;         % meters to feet
         FT_TO_M = 0.3048;          % feet to meters
         M_TO_INCH = 39.37;         % meters to inches
         IN_TO_M = 0.0254;          % inches to meters
         LBF_TO_KG = 0.45359237;    % pounds-force to kilogram
         IN3_TO_M3 = 0.00001639;    % inches cubed to meters cubed

         % Physical Constants
         MU = 3.986004418e14;       % gravitional constant of earth [m^3/s^-2]

         % properties of gases:
         R_U = 8.31446261815324;    % universal gas constant [J·K/mol]

         % Gas properties for air:
         MW_AIR = 28.9647/1000;     % weight of air [kg/mol]
         R_AIR = constant.R_U / constant.MW_AIR;
         GAMMA_AIR = 1.4;

         % Propellant Properties
         LOX_DEN = 0.03922015;      % Density of lox [kg/m^3]
         LNG_DEN = 0.01450439;      % Density of lng [kg/m^3]
         %ETH_DEN =                 % Density of ethanol [kg/m^3]

    end
 end