%RUNFORCDEMO Run a complete FORC analysis against the placeholder solver.
%
%   Set the variables below, then run this script (F5, or type
%   runFORCDemo at the command line). It sweeps the field from Hmax down
%   to Hmin to build the major descending branch, measures a family of
%   first-order reversal curves off it, computes the FORC distribution,
%   and writes everything to OutputFile.
%
%   Every magnetization value comes from simulateHysteresisCurve.m --
%   that is the one file meant to be replaced with a real solver. Nothing
%   below needs to change when it is.
%
%   Once this script has written OutputFile, see plotFORCFamily.m,
%   plotFORCDistribution.m, plotFORCSmoothingMaps.m and
%   plotFORCSmoothingInfluence.m to visualise the results.

%% Parameters -- edit these

Hmax = 1.0;             % top of the field sweep
Hmin = -1.0;            % bottom of the field sweep
NumPoints = 201;        % field steps on the major descending branch

Stride = 5;             % every Stride-th point becomes a reversal field
MeasurementStride = []; % [] keeps the lattice square (= Stride)
HrMin = -Inf;           % restrict reversal fields to this window ...
HrMax = Inf;            % ... (default: the whole sweep)

SmoothingFactor = 3;    % half-width, in grid points, of the FORC distribution fit

OutputFile = fullfile(fileparts(mfilename("fullpath")),"forc_demo_output.json");

%% Run the study

sim = FORCSimulation(Hmax,Hmin, ...
    "NumPoints",NumPoints, ...
    "Stride",Stride, ...
    "MeasurementStride",MeasurementStride, ...
    "HrMin",HrMin, ...
    "HrMax",HrMax, ...
    "SolveFcn",@simulateHysteresisCurve);

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
    "smoothingFactor",SmoothingFactor);
results.majorBranch = struct("H",sim.HMajor,"M",sim.MMajor);
results.diagram = struct("Hr",sim.Hr,"H",sim.H,"M",sim.M);
results.distribution = struct("Hr",dist.Hr,"H",dist.H,"rho",dist.rho, ...
    "SmoothingFactor",dist.SmoothingFactor);

fileID = fopen(OutputFile,"w");
if fileID < 0
    error("runFORCDemo:CannotWrite","Could not open %s for writing.",OutputFile);
end
closeFile = onCleanup(@() fclose(fileID));
fprintf(fileID,"%s",jsonencode(results));

fprintf("Wrote %s\n",OutputFile);
