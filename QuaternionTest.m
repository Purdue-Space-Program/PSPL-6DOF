classdef QuaternionTest < matlab.unittest.TestCase
    % very simpel unit test for quaternions to test CI/CD

    methods (TestClassSetup)
        % Shared setup for the entire test class
    end

    methods (TestMethodSetup)
        % Setup for each test
    end

    methods (Test)
        % Test methods

        function eulvQuat(testCase)
            angle = [.1,0,.1];

            quat = eul2quat(angle, 'ZYX'); % Convert Euler angles to quaternion
            testCase.verifyEqual(quat, [0.9975, 0.049917, 0.0025, 0.049917], 'AbsTol', 1e-5);
        end
    end
end