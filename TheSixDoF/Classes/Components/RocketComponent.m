classdef (Abstract) RocketComponent
    properties
        Name (1,1) string
        Mass (1,1) double
        Position (1,3) double               % CoM or start position [m] (x,y,z)
        Material (1,1) string = "None"
        ComponentLength (1,1) double = 0    % Axial extent [m]; 0 = point mass
        MassDistribution double = []        % Nx2 [norm_pos (0–1), weight]; empty = uniform over ComponentLength
    end

    methods
        function com = getDistributedCoM(obj)
            % Returns [x,y,z] CoM of this component in rocket frame.
            % Position is interpreted as the fore (start) end of the component
            % when ComponentLength > 0, or as the CoM when ComponentLength = 0.
            if obj.ComponentLength == 0 || obj.Mass == 0
                com = obj.Position;
                return
            end

            dist = obj.MassDistribution;
            if isempty(dist)
                s_com = 0.5;    % Uniform: CoM at midpoint
            else
                s   = dist(:,1);
                w   = dist(:,2) / sum(dist(:,2));
                s_com = dot(s, w);
            end

            com    = obj.Position;
            com(1) = obj.Position(1) + s_com * obj.ComponentLength;
        end

        function [Ixx, Iyy, Izz] = getOwnMoI(obj, radius)
            % Returns component's own MoI [Ixx, Iyy, Izz] about its CoM.
            % radius: effective outer radius used for Ixx [m]
            m = obj.Mass;
            L = obj.ComponentLength;

            if m == 0 || L == 0
                Ixx = 0; Iyy = 0; Izz = 0;
                return
            end

            dist = obj.MassDistribution;
            if isempty(dist)
                % Uniform distribution: closed-form thin rod + annulus
                Ixx = 0.5 * m * radius^2;
                Iyy = (1/12) * m * L^2;
                Izz = Iyy;
            else
                % Numerical integration over user-supplied distribution
                s     = dist(:,1);
                w     = dist(:,2) / sum(dist(:,2));
                s_com = dot(s, w);
                ds    = s - s_com;              % Offset from component CoM [normalized]

                Ixx = 0.5 * m * radius^2;
                Iyy = m * dot(w, (ds * L).^2);
                Izz = Iyy;
            end
        end
    end
end
