%RUNFORCDEMOTWOHYSTERONS Run a FORC study of two independent ideal hysterons.
%
%   Set the variables below, then run this script (F5, or type
%   runFORCDemoTwoHysterons at the command line). Plugs in
%   simulateIdealHysteronPopulation.m with two switching fields --
%   Hsw(1) is the same value the single-hysteron demo uses, Hsw(2) is
%   deliberately different -- so the resulting FORC distribution should
%   show TWO separate points at (Hc,Hu) = (Hsw(1),0) and (Hsw(2),0),
%   with nothing connecting them: these particles don't interact, they
%   just happen to be measured together. That clean separation is the
%   whole point of this step -- it's what a later interacting version
%   would change.
%
%   As with the single-hysteron demo, both Hsw values must land exactly
%   on the reversal-field and measurement lattices (see the "grid
%   alignment" note in README.md); the parameters below are chosen so
%   they do.
%
%   Once this script has written OutputFile, call plotFORCFamily or
%   plotFORCDistribution against it (pass OutputPrefix to keep its
%   figures separate from the other demos', and MarkerHcHu with both
%   theoretical points), or run animateFORCIdealHysteron to see both
%   points build up together.

%% Parameters -- edit these

Hmax = 1.0;
Hmin = -1.0;
NumPoints = 201;         % dH = 0.01

Stride = 5;              % reversal-field / measurement spacing = 0.05
MeasurementStride = [];
HrMin = -Inf;
HrMax = Inf;

Hsw = [0.3, 0.6];        % switching field of each particle; both must be
                          % exact multiples of Stride*(Hmax-Hmin)/(NumPoints-1) = 0.05

SmoothingFactor = 2;

OutputFile = fullfile(fileparts(mfilename("fullpath")),"forc_demo_two_hysterons_output.json");

%% Run the study

solveFcn = @(H,Hmax,Hmin,state) simulateIdealHysteronPopulation(H,Hmax,Hmin,state, ...
    "Hsw",Hsw,"Tolerance",1e-9);

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
    "Hsw",Hsw);
results.theoretical = struct("Hc",Hsw(:),"Hu",zeros(numel(Hsw),1));
results.majorBranch = struct("H",sim.HMajor,"M",sim.MMajor);
results.diagram = struct("Hr",sim.Hr,"H",sim.H,"M",sim.M);
results.distribution = struct("Hr",dist.Hr,"H",dist.H,"rho",dist.rho, ...
    "SmoothingFactor",dist.SmoothingFactor);

fileID = fopen(OutputFile,"w");
if fileID < 0
    error("runFORCDemoTwoHysterons:CannotWrite","Could not open %s for writing.",OutputFile);
end
closeFile = onCleanup(@() fclose(fileID));
fprintf(fileID,"%s",jsonencode(results));

fprintf("Wrote %s\n",OutputFile);
