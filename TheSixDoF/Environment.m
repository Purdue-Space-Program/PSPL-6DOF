classdef Environment
    % The environment class models the environment system for the rocket,
    % including the sea level height of the location, local gravity, 

    properties
        date (1,1) datetime = datetime("today");
        Elevation (1,1) double = 614;
        latlong (1,2) double = [35.3462082, -117.8166643]
        geocentricRadius (1,1) double = 6.371077849286893e6;
    end

    methods
        % this shit not working rn
        % function elev = getElevation(lat,long)
        %     loc = txsite("Latitude", lat,"Longitude",long);
        %     elev = elevation(loc);
        % end
    end
end