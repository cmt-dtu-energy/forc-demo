classdef testForcDistribution < matlab.unittest.TestCase
    properties (TestParameter)
        smoothingFactor = {1, 2, 3};
        gridSizeNumbers = {9, 11, 13};
    end

    methods (Test)
        function testBilinearSurfaceGivesKnownDerivative(testCase, smoothingFactor)
            n = 17;
            x = linspace(-1, 1, n);
            y = linspace(-1, 1, n);
            [Y, X] = meshgrid(y, x);
            Z = X .* Y;
            d = forcDistribution(x, y, Z, smoothingFactor);

            testCase.verifySize(d.rho, size(Z));
            testCase.verifyEqual(d.('Hr'), x);
            testCase.verifyEqual(d.('H'),  y);

            interior = ~isnan(d.rho);
            expected = repmat(-0.5, nnz(interior), 1);
            testCase.verifyEqual(d.rho(interior), expected, 'AbsTol', 1e-12);
        end

        function testEdgeWindowLeavesEdgesUndefined(testCase, smoothingFactor, gridSizeNumbers)
            x = linspace(-1, 1, gridSizeNumbers);
            y = linspace(-1, 1, gridSizeNumbers);
            [Y, X] = meshgrid(y, x);
            Z = X .* Y;
            d = forcDistribution(x, y, Z, smoothingFactor);

            expectedUndefined = false(gridSizeNumbers, gridSizeNumbers);
            expectedUndefined([1:smoothingFactor, gridSizeNumbers-smoothingFactor+1:gridSizeNumbers], :) = true;
            expectedUndefined(:, [1:smoothingFactor, gridSizeNumbers-smoothingFactor+1:gridSizeNumbers]) = true;
            testCase.verifyTrue(all(isnan(d.rho(expectedUndefined))));
            testCase.verifyTrue(all(~isnan(d.rho(~expectedUndefined))));
        end

        function testQuadraticInFirstCoordinateHasZeroMixedDerivative(testCase, smoothingFactor)
            n = 17;
            x = linspace(-1, 1, n);
            y = linspace(-1, 1, n);
            [~, X] = meshgrid(y, x);
            Z = X.^2;
            d = forcDistribution(x, y, Z, smoothingFactor);

            interior = ~isnan(d.rho);
            expected = zeros(nnz(interior), 1);
            testCase.verifyEqual(d.rho(interior), expected, 'AbsTol', 1e-12);
        end
    end
end
