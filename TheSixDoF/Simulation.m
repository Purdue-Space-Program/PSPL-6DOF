classdef Simulation
    % The simulation class sets all of the settings for the simulation of
    % the rocket. Various parameters can be used to change the length and
    % fidelity of the simulation

    properties
        %EndCondition (1,1) string {mustBeMember(EndCondition,["apogee","full","burnout"])} = "apogee"
        EndCondition (1,1) string = 'apogee'

        Output (1,1) logical = 1
        Timestep (1,1) double = 0.1
        %Fidelity (1,1) string {mustBeMember(Fidelity,["low","medium","high"])} = "medium"
        Fidelity (1,1) string = "medium"

    end

    methods

        function sim = Simulation(endCondition, timestep, fidelity, output)

            if (nargin == 4)
                endCondition = convertCharsToStrings(endCondition);
                fidelity = convertCharsToStrings(fidelity);

                if isstring(endCondition)
                    sim.EndCondition = endCondition;
                else
                    error("The end condition must be text")
                end
                sim.EndCondition = endCondition;
                sim.Timestep = timestep;
                sim.Fidelity = fidelity;
                sim.Output = output;
            end

        end
    end
        methods (Access = private)

        end

    end
