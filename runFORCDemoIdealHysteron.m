%RUNFORCDEMOIDEALHYSTERON Run a FORC study of a single ideal hysteron.
%
%   Set the variables below, then run this script (F5, or type
%   runFORCDemoIdealHysteron at the command line). Unlike runFORCDemo.m,
%   which exercises the placeholder smoothed solver, this run plugs in
%   simulateIdealHysteron.m -- a particle that switches instantaneously
%   at +-Hsw -- so the resulting FORC distribution should concentrate
%   near a single point at (Hc,Hu) = (Hsw,0) rather than a smoothed
%   ridge. That only comes out cleanly if Hsw lands exactly on both the
%   reversal-field lattice and the measurement lattice; the parameters
%   below are chosen so it does (see the "grid alignment" note in
%   README.md).
%
%   Even with perfect alignment, don't expect a single-cell spike:
%   because the switching corner then sits exactly on a grid vertex
%   shared by four cells, Pike's fit (see forcDistribution.m) responds
%   equally at all four, giving a small symmetric plateau centered half
%   a coercivity grid-step away from the theoretical point rather than
%   a single value exactly on it. This is intrinsic to a centered fit
%   resolving an on-lattice singularity, not a bug -- see README.md.
%
%   Once this script has written OutputFile, call plotFORCFamily,
%   plotFORCDistribution, plotFORCSmoothingMaps or
%   plotFORCSmoothingInfluence against it (pass OutputPrefix to keep
%   its figures separate from runFORCDemo.m's), or run
%   animateFORCIdealHysteron to see the point build up.

%% Parameters -- edit these

Hmax = 1.0;             % top of the field sweep
Hmin = -1.0;            % bottom of the field sweep
NumPoints = 201;        % field steps on the major descending branch, dH = 0.01

Stride = 5;             % every Stride-th point becomes a reversal field
MeasurementStride = []; % [] keeps the lattice square (= Stride)
HrMin = -Inf;           % restrict reversal fields to this window ...
HrMax = Inf;            % ... (default: the whole sweep)

Hsw = 0.3;              % switching field magnitude of the ideal hysteron.
                         % With the grid above, the reversal-field and
                         % measurement spacing is Stride*(Hmax-Hmin)/(NumPoints-1)
                         % = 5*0.01 = 0.05, and Hsw = 6*0.05 lands exactly
                         % on that lattice -- change Hsw only together
                         % with the grid parameters, keeping it an exact
                         % multiple of that spacing.

SmoothingFactor = 2;    % half-width, in grid points, of the FORC distribution fit;
                         % kept tighter than runFORCDemo.m's 3 so the
                         % near-delta stays as sharp as the fit allows

OutputFile = fullfile(fileparts(mfilename("fullpath")),"forc_demo_ideal_hysteron_output.json");

%% Run the study

solveFcn = @(H,Hmax,Hmin,state) simulateIdealHysteron(H,Hmax,Hmin,state, ...
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
results.theoretical = struct("Hc",Hsw,"Hu",0);
results.majorBranch = struct("H",sim.HMajor,"M",sim.MMajor);
results.diagram = struct("Hr",sim.Hr,"H",sim.H,"M",sim.M);
results.distribution = struct("Hr",dist.Hr,"H",dist.H,"rho",dist.rho, ...
    "SmoothingFactor",dist.SmoothingFactor);

fileID = fopen(OutputFile,"w");
if fileID < 0
    error("runFORCDemoIdealHysteron:CannotWrite","Could not open %s for writing.",OutputFile);
end
closeFile = onCleanup(@() fclose(fileID));
fprintf(fileID,"%s",jsonencode(results));

fprintf("Wrote %s\n",OutputFile);
