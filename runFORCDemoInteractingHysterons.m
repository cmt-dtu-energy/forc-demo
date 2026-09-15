%RUNFORCDEMOINTERACTINGHYSTERONS Run a FORC study of two dipolar-coupled ideal hysterons.
%
%   Set the variables below, then run this script. Plugs in
%   simulateInteractingHysterons.m: the same two switching fields as
%   runFORCDemoTwoHysterons.m (Hsw = [0.3, 0.6]), but now coupled by a
%   dipolar interaction field, CouplingField. Positive CouplingField is
%   the "chain" (head-to-tail) geometry -- favours parallel alignment,
%   so once one particle switches it helps drag the other along. See
%   simulateInteractingHysterons.m for the sign convention and the
%   dipole-field derivation behind it.
%
%   Every Hsw(i) +- CouplingField needs to land on the reversal-field/
%   measurement lattice for the same reason a bare Hsw does in the
%   single/two-hysteron demos (see the "grid alignment" note in
%   README.md) -- CouplingField below is an exact multiple of the grid
%   spacing so every shifted threshold stays on-grid too.
%
%   Once this script has written OutputFile, call plotFORCFamily,
%   plotFORCDistribution, plotFORCSmoothingMaps or
%   plotFORCSmoothingInfluence against it, or run
%   animateFORCIdealHysteron to see it build up. Compare directly
%   against forc_demo_two_hysterons_output.json (same Hsw, no
%   coupling) to see what the interaction actually changed -- see
%   README.md for what was observed with the parameters below.

%% Parameters -- edit these

Hmax = 1.0;
Hmin = -1.0;
NumPoints = 201;         % dH = 0.01

Stride = 5;              % reversal-field / measurement spacing = 0.05
MeasurementStride = [];
HrMin = -Inf;
HrMax = Inf;

Hsw = [0.3, 0.6];             % each particle's own (uncoupled) switching field
CouplingField = 0.1;          % dipolar coupling; positive = chain/ferromagnetic-like.
                               % An exact multiple of the grid spacing (0.05) so every
                               % Hsw(i) +- CouplingField stays on-grid; below the
                               % coupling, (Hsw(2)-Hsw(1))/2 = 0.15, that would make
                               % the two ascending thresholds coincide (see
                               % simulateInteractingHysterons.m).

SmoothingFactor = 2;

OutputFile = fullfile(fileparts(mfilename("fullpath")),"forc_demo_interacting_hysterons_output.json");

%% Run the study

solveFcn = @(H,Hmax,Hmin,state) simulateInteractingHysterons(H,Hmax,Hmin,state, ...
    "Hsw",Hsw,"CouplingField",CouplingField,"Tolerance",1e-9);

sim = FORCSimulation(Hmax,Hmin, ...
    "NumPoints",NumPoints, ...
    "Stride",Stride, ...
    "MeasurementStride",MeasurementStride, ...
    "HrMin",HrMin, ...
    "HrMax",HrMax, ...
    "SolveFcn",solveFcn);

sim.run();
dist = sim.distribution(SmoothingFactor);

%% Write the results

results.parameters = struct( ...
    "Hmax",Hmax, ...
    "Hmin",Hmin, ...
    "numPoints",NumPoints, ...
    "stride",Stride, ...
    "HrMin",HrMin, ...
    "HrMax",HrMax, ...
    "smoothingFactor",SmoothingFactor, ...
    "Hsw",Hsw, ...
    "CouplingField",CouplingField);
% No "theoretical" field: unlike the non-interacting demos, the exact
% peak location here depends on coupled switching history and isn't a
% simple closed form -- see README.md for what was actually observed.
% isolatedReference records where the two points sit with no coupling
% at all (CouplingField=0, i.e. the forc_demo_two_hysterons_output.json
% case), purely as a visual comparison.
results.isolatedReference = struct("Hc",Hsw(:),"Hu",zeros(numel(Hsw),1));
results.majorBranch = struct("H",sim.HMajor,"M",sim.MMajor);
results.diagram = struct("Hr",sim.Hr,"H",sim.H,"M",sim.M);
results.distribution = struct("Hr",dist.Hr,"H",dist.H,"rho",dist.rho, ...
    "SmoothingFactor",dist.SmoothingFactor);

fileID = fopen(OutputFile,"w");
if fileID < 0
    error("runFORCDemoInteractingHysterons:CannotWrite","Could not open %s for writing.",OutputFile);
end
closeFile = onCleanup(@() fclose(fileID));
fprintf(fileID,"%s",jsonencode(results));

fprintf("Wrote %s\n",OutputFile);
