classdef Simulation
    % The simulation class sets all of the settings for the simulation of
    % the rocket. Various parameters can be used to change the length and
    % fidelity of the simulation

    properties
        EndCondition (1,1) string {mustBeMember(EndCondition,["apogee","full","burnout"])} = "apogee"
        Output (1,1) logical = 1
        Timestep (1,1) double = 0.1
        Fidelity (1,1) string {mustBeMember(Fidelity,["low","medium","high"])} = "medium"
    end

    properties (SetAccess = private)
        relTol (1,1) double = 1e-6;
        absTol (1,1) double = 1e-6;
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

                % set the fidelity options for the sim
                sim = SetIntegrationParams(sim);
            end

        end
    end
    methods (Access = private)
        function sim = SetIntegrationParams(sim)
            fidelity = sim.Fidelity;
            switch fidelity
                case "low"
                    sim.relTol = 1e-2;
                    sim.absTol = 1e-2;
                case "medium"
                    sim.relTol = 1e-3;
                    sim.absTol = 1e-6;
                case "high"
                    sim.relTol = 1e-4;
                    sim.absTol = 1e-7;
            end
        end
    end
end
