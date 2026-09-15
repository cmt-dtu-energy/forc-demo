classdef testForcDistribution < matlab.unittest.TestCase
    methods (Test)
        function bilinearSurfaceGivesKnownDerivative(testCase)
            smoothingWidth = 2;
            n = 17;
            x = linspace(-1, 1, n);
            y = linspace(-1, 1, n);
            [Y, X] = meshgrid(y, x);
            Z = X .* Y;
            d = forcDistribution(x, y, Z, smoothingWidth);
            testCase.verifySize(d.rho, size(Z));
            testCase.verifyEqual(d.('Hr'), x);
            testCase.verifyEqual(d.('H'),  y);
            interior = ~isnan(d.rho);
            expected = repmat(-0.5, nnz(interior), 1);
            testCase.verifyEqual(d.rho(interior), expected, 'AbsTol', 1e-12);
        end

        function edgeWindowLeavesEdgesUndefined(testCase)
            smoothingWidth = 1;
            n = 9;
            x = linspace(-1, 1, n);
            y = linspace(-1, 1, n);
            [Y, X] = meshgrid(y, x);
            Z = X .* Y;
            d = forcDistribution(x, y, Z, smoothingWidth);
            expectedUndefined = false(n, n);
            expectedUndefined([1:smoothingWidth, n-smoothingWidth+1:n], :) = true;
            expectedUndefined(:, [1:smoothingWidth, n-smoothingWidth+1:n]) = true;
            testCase.verifyTrue(all(isnan(d.rho(expectedUndefined))));
            testCase.verifyTrue(all(~isnan(d.rho(~expectedUndefined))));
        end

        function quadraticInFirstCoordinateHasZeroMixedDerivative(testCase)
            smoothingWidth = 2;
            n = 17;
            x = linspace(-1, 1, n);
            y = linspace(-1, 1, n);
            [~, X] = meshgrid(y, x);
            Z = X.^2;
            d = forcDistribution(x, y, Z, smoothingWidth);
            interior = ~isnan(d.rho);
            expected = zeros(nnz(interior), 1);
            testCase.verifyEqual(d.rho(interior), expected, 'AbsTol', 1e-12);
        end

        function nonUniformGridRecoversKnownDerivative(testCase)
            smoothingWidth = 2;
            nXDir = 21;
            nYDir = 19;
            x = linspace(-1, 1, nXDir).^3;
            y = logspace(-2, 0, nYDir);
            [Y, X] = meshgrid(y, x);
            Z = X .* Y;
            d = forcDistribution(x, y, Z, smoothingWidth);
            interior = ~isnan(d.rho);
            expected = repmat(-0.5, nnz(interior), 1);
            testCase.verifyEqual(d.rho(interior), expected, 'AbsTol', 1e-10);
        end
    end
end
