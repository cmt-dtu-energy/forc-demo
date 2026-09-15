classdef testForcDistribution < matlab.unittest.TestCase
    %TESTFORCDISTRIBUTION Unit tests for the shared FORC distribution fit.

    methods (Test)
        function bilinearSurfaceGivesKnownDerivative(testCase)
            % For M(rev, meas) = rev * meas, d2M/(drev dmeas) = 1,
            % so rho should be -0.5 at every interior grid point.

            smoothingWidth = 2;
            n = 17;
            revFields = linspace(-1, 1, n);
            measFields = linspace(-1, 1, n);
            [measGrid, revGrid] = meshgrid(measFields, revFields);
            fieldMatrix = revGrid .* measGrid;

            d = forcDistribution(revFields, measFields, fieldMatrix, smoothingWidth);

            testCase.verifySize(d.rho, size(fieldMatrix));
            testCase.verifyEqual(d.('Hr'), revFields);
            testCase.verifyEqual(d.('H'),  measFields);

            interior = ~isnan(d.rho);
            testCase.verifyTrue(any(interior(:)));

            expected = repmat(-0.5, nnz(interior), 1);
            testCase.verifyEqual(d.rho(interior), expected, 'AbsTol', 1e-12);
        end

        function smoothingWindowLeavesEdgesUndefined(testCase)
            % Points within SF rows/columns of the border cannot fit a
            % full (2*SF+1)^2 neighbourhood, so they must stay undefined.

            smoothingWidth = 1;
            n = 9;
            revFields = linspace(-1, 1, n);
            measFields = linspace(-1, 1, n);
            [measGrid, revGrid] = meshgrid(measFields, revFields);
            fieldMatrix = revGrid .* measGrid;

            d = forcDistribution(revFields, measFields, fieldMatrix, smoothingWidth);

            expectedUndefined = false(n, n);
            expectedUndefined([1:smoothingWidth, n-smoothingWidth+1:n], :) = true;
            expectedUndefined(:, [1:smoothingWidth, n-smoothingWidth+1:n]) = true;

            testCase.verifyTrue(all(isnan(d.rho(expectedUndefined))));
            testCase.verifyTrue(all(~isnan(d.rho(~expectedUndefined))));
        end
    end
end
