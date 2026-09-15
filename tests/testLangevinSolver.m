classdef testLangevinSolver < matlab.unittest.TestCase
    methods (Test)
        function zeroFieldHasZeroMagnetization(testCase)
            [M,state] = simulateLangevin(0,1,-1,struct("marker",7));
            testCase.verifyEqual(M,0,"AbsTol",1e-14);
            testCase.verifyEqual(state,struct("marker",7));
        end

        function magnetizationIsOddAndMonotone(testCase)
            fields = linspace(-1,1,101);
            values = arrayfun(@(H) simulateLangevin(H,1,-1,[]),fields);
            testCase.verifyEqual(values,-fliplr(values),"AbsTol",1e-14);
            testCase.verifyGreaterThanOrEqual(diff(values),-1e-14);
            testCase.verifyGreaterThanOrEqual(values,-1);
            testCase.verifyLessThanOrEqual(values,1);
        end

        function parametersControlScaleAndOffset(testCase)
            [positive,~] = simulateLangevin(1,2,-2,[],"Ms",2,"Alpha",3);
            [negative,~] = simulateLangevin(-1,2,-2,[],"Ms",2,"Alpha",3);
            [shifted,~] = simulateLangevin(0,2,-2,[],"Ms",2,"Alpha",3,"BiasField",1);
            testCase.verifyEqual(positive,-negative,"AbsTol",1e-14);
            testCase.verifyGreaterThan(positive,1.3);
            testCase.verifyEqual(shifted,positive,"AbsTol",1e-14);
        end

        function largeFieldsApproachSaturation(testCase)
            [positive,~] = simulateLangevin(1e6,1e6 + 1,-1e6 - 1,[]);
            [negative,~] = simulateLangevin(-1e6,1e6 + 1,-1e6 - 1,[]);
            testCase.verifyEqual(positive,1,"AbsTol",1e-6);
            testCase.verifyEqual(negative,-1,"AbsTol",1e-6);
        end

        function smallArgumentsUseStableLimit(testCase)
            [M,~] = simulateLangevin(1e-8,1,-1,[],"Alpha",2);
            testCase.verifyEqual(M,2e-8/3,"AbsTol",1e-15);
        end

        function invalidParametersAreRejected(testCase)
            testCase.verifyError(@() simulateLangevin(0,1,-1,[],"Ms",0), ...
                "MATLAB:validators:mustBePositive");
            testCase.verifyError(@() simulateLangevin(0,1,-1,[],"Alpha",0), ...
                "MATLAB:validators:mustBePositive");
        end

        function invalidWindowIsRejected(testCase)
            testCase.verifyError(@() simulateLangevin(0,-1,1,[]), ...
                "simulateLangevin:BadWindow");
        end

        function solverIntegratesWithForcSimulation(testCase)
            solveFcn = @(H,Hmax,Hmin,state) ...
                simulateLangevin(H,Hmax,Hmin,state,"Alpha",4);
            study = FORCSimulation(1,-1,"NumPoints",41,"Stride",4, ...
                "SolveFcn",solveFcn,"Verbose",false);
            study.run();
            expected = arrayfun(@(H) simulateLangevin(H,1,-1,[],"Alpha",4),study.H);
            testCase.verifyFalse(isempty(study.M));
            for i = 1:size(study.M,1)
                columns = ~isnan(study.M(i,:));
                testCase.verifyEqual(study.M(i,columns),expected(columns),"AbsTol",1e-12);
            end
        end
    end
end
