function [M,state] = simulateHysteresisCurve(H,Hmax,Hmin,state)
%SIMULATEHYSTERESISCURVE Placeholder single-step magnetization evaluator.
%
%   [M,STATE] = SIMULATEHYSTERESISCURVE(H,HMAX,HMIN,STATE) returns the
%   magnetization M at applied field H, given the saturation field
%   window [HMIN, HMAX] and the solver's previous STATE.
%
%   STATE is opaque to every other file in this bundle -- it is just
%   threaded through repeated calls the same way a real physical solver
%   would carry forward its own memory of the sample (a magnetization
%   direction, a set of internal variables, whatever the model needs).
%   It is also exactly what gets snapshotted at every field step of the
%   descending curve and handed back unchanged at each FORC reversal
%   point, so FORCEngine never has to know what is inside it.
%
%   Call SIMULATEHYSTERESISCURVE(H,HMAX,HMIN,[]) to start a fresh curve
%   at saturation; an empty STATE is treated as "no history yet".
%
%   THIS IS A PLACEHOLDER. Replace the body below with a real
%   single-point magnetization solver -- nothing else in this bundle
%   (FORCFieldPlan, FORCEngine, FORCDiagram, FORCDistribution, or the
%   plotting scripts) knows or cares what happens inside this function;
%   they all work purely off the H/M arrays it produces.
%
%   The stub implemented here is a single hysteron with coercive field
%   Hc: it lags the saturating tanh curve M = Ms*tanh(H/Hk) by Hc in
%   whichever direction the field last moved, which is just enough
%   memory to produce a family of nested minor loops for the rest of the
%   pipeline to be exercised against.

arguments (Input)
    H (1,1) double {mustBeFinite, mustBeReal}
    Hmax (1,1) double {mustBeFinite, mustBeReal}
    Hmin (1,1) double {mustBeFinite, mustBeReal}
    state
end

arguments (Output)
    M (1,1) double
    state
end

Ms = 1;
Hk = 0.4*(Hmax - Hmin);
Hc = 0.15*(Hmax - Hmin);

if isempty(state)
    state = struct("Hprev",Hmax,"direction",-1);
end

direction = sign(H - state.Hprev);
if direction == 0
    direction = state.direction;
end

Heff = H - direction*Hc;
M = Ms*tanh(Heff/Hk);

state.Hprev = H;
state.direction = direction;

end
