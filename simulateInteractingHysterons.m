function [M,state] = simulateInteractingHysterons(H,Hmax,Hmin,state,options)
%SIMULATEINTERACTINGHYSTERONS Two ideal hysterons coupled by a dipolar interaction field.
%
%   [M,STATE] = SIMULATEINTERACTINGHYSTERONS(H,HMAX,HMIN,STATE) evolves
%   two ideal rectangular hysterons (see simulateIdealHysteron.m) that
%   are no longer independent: each one's effective field is the
%   applied field H plus a term proportional to the OTHER particle's
%   current magnetization,
%
%       Heff_i = H + CouplingField * M_j      (j ~= i)
%
%   which is the standard mean-field simplification of a dipolar
%   interaction. For two point dipoles a fixed distance r apart, the
%   field one exerts on the other along their shared axis is
%   B = (mu0/4pi)*[3(m.rhat)rhat - m]/r^3: for a "chain" (head-to-tail)
%   geometry this works out to a field in the SAME direction as the
%   source dipole (favours parallel alignment -- ferromagnetic-like);
%   for a "side-by-side" geometry it flips sign (favours antiparallel
%   alignment -- antiferromagnetic-like). CouplingField folds that
%   whole geometric dependence into one signed number: positive for a
%   chain-like arrangement, negative for side-by-side, magnitude set by
%   the real mu0*m/(4pi*r^3) at whatever separation is being modelled.
%   This repo never tracks real 3D positions (every other solver here
%   is a scalar hysteron too), so a fixed separation/orientation is
%   assumed rather than simulated explicitly.
%
%   Each particle's threshold-crossing check uses the OTHER's state as
%   of the END of the previous field step, not a value updated within
%   the same step (a staggered/Jacobi update) -- this keeps the result
%   deterministic and free of any particle-ordering choice, and is
%   exact in the limit of a field grid fine enough that both particles
%   cannot cross their thresholds at truly the same instant (this
%   repo's usual dH is far finer than any Hsw or CouplingField used
%   here).
%
%   [...] = SIMULATEINTERACTINGHYSTERONS(...,"Hsw",HSW) sets each
%   particle's own (uncoupled) switching-field magnitude, a 1x2 vector
%   (default [0.3 0.3]; pass distinct values to tell them apart).
%
%   [...] = SIMULATEINTERACTINGHYSTERONS(...,"CouplingField",J) sets
%   the interaction strength and sign (default 0.1). J=0 recovers
%   simulateIdealHysteronPopulation.m's independent-particle behaviour
%   exactly.
%
%   See also simulateIdealHysteron, simulateIdealHysteronPopulation,
%   runFORCDemoInteractingHysterons.

arguments (Input)
    H (1,1) double {mustBeFinite, mustBeReal}
    Hmax (1,1) double {mustBeFinite, mustBeReal}
    Hmin (1,1) double {mustBeFinite, mustBeReal}
    state
    options.Hsw (1,2) double {mustBePositive} = [0.3 0.3]
    options.CouplingField (1,1) double {mustBeFinite, mustBeReal} = 0.1
    options.Tolerance (1,1) double {mustBeNonnegative} = 1e-9
end

arguments (Output)
    M (1,1) double
    state
end

if isempty(state)
    state = {[],[]};
end

M1prev = currentM_(state{1});
M2prev = currentM_(state{2});

Heff1 = H + options.CouplingField*M2prev;
Heff2 = H + options.CouplingField*M1prev;

[M1,state{1}] = simulateIdealHysteron(Heff1,Hmax,Hmin,state{1}, ...
    "Hsw",options.Hsw(1),"Tolerance",options.Tolerance);
[M2,state{2}] = simulateIdealHysteron(Heff2,Hmax,Hmin,state{2}, ...
    "Hsw",options.Hsw(2),"Tolerance",options.Tolerance);

M = mean([M1,M2]);

end


function m = currentM_(s)
%CURRENTM_ The magnetization a hysteron's state encodes, treating empty as fresh (+1).
if isempty(s)
    m = 1;
else
    m = s;
end
end
