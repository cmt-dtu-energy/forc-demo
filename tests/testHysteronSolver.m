classdef testHysteronSolver < matlab.unittest.TestCase
    methods (Test)
        function solverRunsFullSimulation(testCase)
            upperLimit = 1.0;
            lowerLimit = -1.0;
            study = FORCSimulation(upperLimit, lowerLimit, ...
                'NumPoints', 101, ...
                'Stride', 5, ...
                'Verbose', false);
            study.run();

            testCase.verifyFalse(isempty(study.M));
            testCase.verifyTrue(any(~isnan(study.M(:))));
            testCase.verifyEqual(study.HMajor(1), upperLimit, 'AbsTol', 1e-12);
            testCase.verifyEqual(study.HMajor(end), lowerLimit, 'AbsTol', 1e-12);
        end

        function solverProducesBoundedMagnetization(testCase)
            upperLimit = 1.0;
            lowerLimit = -1.0;
            study = FORCSimulation(upperLimit, lowerLimit, ...
                'NumPoints', 101, ...
                'Stride', 5, ...
                'Verbose', false);
            study.run();

            measuredValues = study.M(~isnan(study.M));
            testCase.verifyGreaterThanOrEqual(measuredValues, -1 - 1e-12);
            testCase.verifyLessThanOrEqual(measuredValues, 1 + 1e-12);
        end

        function strideAlignsReversals(testCase)
            upperLimit = 1.0;
            lowerLimit = -1.0;
            strideValue = 5;
            study = FORCSimulation(upperLimit, lowerLimit, ...
                'NumPoints', 101, ...
                'Stride', strideValue, ...
                'Verbose', false);
            study.run();

            expectedReversals = study.HMajor(1 + strideValue : strideValue : end);
            testCase.verifyEqual(study.Hr, expectedReversals, 'AbsTol', 1e-12);
        end
    end
end
