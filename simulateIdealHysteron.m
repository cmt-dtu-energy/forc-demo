function [M,state] = simulateIdealHysteron(H,Hmax,Hmin,state,options)
%SIMULATEIDEALHYSTERON An ideal rectangular (Preisach) hysteron.
%
%   [M,STATE] = SIMULATEIDEALHYSTERON(H,HMAX,HMIN,STATE) returns the
%   magnetization M at applied field H for a single particle that
%   switches instantaneously between -1 and +1: it flips up at
%   H = +Hsw and down at H = -Hsw, and holds its value everywhere else.
%   Hsw is an unsigned magnitude -- the two switching fields +Hsw and
%   -Hsw are both derived from it here, which is what makes this
%   hysteron symmetric (no bias field). This is the textbook Preisach
%   unit: its FORC distribution is a single delta function at
%   (Hc,Hu) = (Hsw,0), rather than a smoothed ridge.
%
%   Call SIMULATEIDEALHYSTERON(H,HMAX,HMIN,[]) to start a fresh curve
%   at saturation; an empty STATE is treated as "starts at Hmax,
%   saturated positive" -- the same convention FORCSimulation's major
%   branch relies on.
%
%   [...] = SIMULATEIDEALHYSTERON(...,"Hsw",HSW) sets the switching
%   field magnitude (default 0.3). For a demo that needs the switching
%   delta to land exactly on the FORC grid, HSW must be an exact
%   multiple of both the reversal-field spacing and the measurement
%   spacing (see runFORCDemoIdealHysteron.m).
%
%   [...] = SIMULATEIDEALHYSTERON(...,"Tolerance",TOL) sets the
%   half-width, in field units, of the crossing test around +-Hsw
%   (default 1e-9). This is not cosmetic: HMajor = linspace(Hmax,Hmin,
%   NumPoints) is not guaranteed to place a grid point at exactly the
%   literal value of Hsw down to the last bit, and without slack here
%   the switching event could land one grid cell away from where
%   Hsw says it should, silently smearing the delta that is the whole
%   point of this demo.
%
%   See also FORCSimulation, simulateHysteresisCurve.

arguments (Input)
    H (1,1) double {mustBeFinite, mustBeReal}
    Hmax (1,1) double {mustBeFinite, mustBeReal}
    Hmin (1,1) double {mustBeFinite, mustBeReal}
    state
    options.Hsw (1,1) double {mustBePositive} = 0.3
    options.Tolerance (1,1) double {mustBeNonnegative} = 1e-9
end

arguments (Output)
    M (1,1) double
    state
end

if options.Hsw >= min(Hmax,-Hmin)
    error("simulateIdealHysteron:BadHsw", ...
          "Hsw (%g) must be below both Hmax (%g) and -Hmin (%g), so that both +Hsw and -Hsw lie strictly inside [Hmin, Hmax].", ...
          options.Hsw,Hmax,-Hmin);
end

if isempty(state)
    state = 1;
end

if H >= options.Hsw - options.Tolerance
    state = 1;
elseif H <= -options.Hsw + options.Tolerance
    state = -1;
end

M = state;

end
