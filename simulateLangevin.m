function [M,state] = simulateLangevin(H,Hmax,Hmin,state,options)
arguments (Input)
    H (1,1) double {mustBeFinite, mustBeReal}
    Hmax (1,1) double {mustBeFinite, mustBeReal}
    Hmin (1,1) double {mustBeFinite, mustBeReal}
    state
    options.Ms (1,1) double {mustBePositive, mustBeFinite} = 1
    options.Alpha (1,1) double {mustBePositive, mustBeFinite} = 4
    options.BiasField (1,1) double {mustBeFinite, mustBeReal} = 0
end
arguments (Output)
    M (1,1) double
    state
end
if Hmax <= Hmin
    error("simulateLangevin:BadWindow", ...
        "Hmax (%g) must be above Hmin (%g).",Hmax,Hmin);
end
x = options.Alpha*(H + options.BiasField);
M = options.Ms*langevin_(x);
end

function value = langevin_(x)
small = abs(x) < 1e-4;
value = zeros(size(x));
value(~small) = 1./tanh(x(~small)) - 1./x(~small);
value(small) = x(small)/3 - x(small).^3/45 + 2*x(small).^5/945;
end
