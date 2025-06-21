classdef Environment
    % The environment class models the environment system for the rocket,
    % including the sea level height of the location, local gravity, and

    properties
        Date (1,1) datetime = datetime("today");
        Elevation (1,1) double = 627.91;
        LatLong (1,2) double = [35.347444074690735, -117.8090720168799]
        geocentricRadius (1,1) double = 6.371077849286893e6;
    end

    methods
        function env = setLaunchSite(env, lat, long, date)
            % setLaunchSite creates a new launch site given an input
            % environment and the position in latitude and longitude. The
            % function automatically generates
            if (nargin == 3)
                env = Env.Environment;
                env.LatLong = [lat,long];
                env.Elevation = getElevation(env);
                env.Date = env.Date;
            elseif (nargin == 4)
                env.Date = date;
            end
        end

        function elev = getElevation(env)
            %env = Environment;
            loc = txsite("Latitude", env.LatLong(1),"Longitude",env.LatLong(2));
            elev = elevation(loc);
        end
    end
end