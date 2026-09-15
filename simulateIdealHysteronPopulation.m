function [M,state] = simulateIdealHysteronPopulation(H,Hmax,Hmin,state,options)
%SIMULATEIDEALHYSTERONPOPULATION A set of independent, non-interacting ideal hysterons.
%
%   [M,STATE] = SIMULATEIDEALHYSTERONPOPULATION(H,HMAX,HMIN,STATE)
%   evolves several ideal rectangular hysterons (see
%   simulateIdealHysteron.m) side by side, each switching only on the
%   applied field H -- never on the others' state, so there is no
%   interaction between them yet. M is their normalized average (each
%   particle contributes +-1/N, so the population still saturates at
%   +-1, just like a single hysteron), which keeps this solver's M
%   directly comparable to simulateIdealHysteron.m's.
%
%   [...] = SIMULATEIDEALHYSTERONPOPULATION(...,"Hsw",HSW) sets the
%   vector of switching-field magnitudes, one per particle (default
%   [0.3 0.3], i.e. two particles, kept as a placeholder default; pass
%   distinct values to tell them apart). A single-element HSW
%   degenerates to exactly simulateIdealHysteron.m's behaviour.
%
%   With no interaction between particles, this population's FORC
%   distribution is the SUM of each particle's own delta -- one clean
%   point at (Hc,Hu) = (Hsw(i),0) per particle, with no ridge or
%   off-axis structure connecting them. That absence of coupling is
%   the point of this demo: it is exactly what a later interacting
%   version would change.
%
%   See also simulateIdealHysteron, runFORCDemoTwoHysterons.

arguments (Input)
    H (1,1) double {mustBeFinite, mustBeReal}
    Hmax (1,1) double {mustBeFinite, mustBeReal}
    Hmin (1,1) double {mustBeFinite, mustBeReal}
    state
    options.Hsw (1,:) double {mustBePositive} = [0.3 0.3]
    options.Tolerance (1,1) double {mustBeNonnegative} = 1e-9
end

arguments (Output)
    M (1,1) double
    state
end

n = numel(options.Hsw);

if isempty(state)
    state = cell(1,n);
end

Mi = zeros(1,n);
for i = 1:n
    [Mi(i),state{i}] = simulateIdealHysteron(H,Hmax,Hmin,state{i}, ...
        "Hsw",options.Hsw(i),"Tolerance",options.Tolerance);
end

M = mean(Mi);

end
