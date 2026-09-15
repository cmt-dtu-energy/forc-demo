classdef FORCSimulation < handle
    %FORCSIMULATION A First-Order Reversal Curve (FORC) study of a sample.
    %
    %   A FORC study repeats one experiment many times: bring the sample
    %   down the major descending branch of its hysteresis loop to some
    %   "reversal field" Hr, then sweep the field back up to saturation
    %   and record the magnetization on the way. Every reversal curve
    %   becomes one row of a family M(Hr,H); its mixed second derivative,
    %
    %       rho(Hr,H) = -1/2 d2M / (dHr dH)
    %
    %   is the FORC distribution. A single hysteron (one particle
    %   switching irreversibly at a single field) contributes a delta
    %   function to rho, so the distribution reads as the density of
    %   switching events: its spread along the coercivity coordinate
    %   Hc = (H-Hr)/2 is the spread of coercivities in the sample, and
    %   its spread along the interaction coordinate Hu = (H+Hr)/2 is the
    %   spread of interaction fields between particles.
    %
    %   Every magnetization value used here comes from one function
    %   handle, SolveFcn -- by default simulateHysteresisCurve.m -- called
    %   as
    %
    %       [M,state] = SolveFcn(H,Hmax,Hmin,state)
    %
    %   with STATE threaded from one call to the next. STATE is opaque to
    %   this class: it is just recorded after every step of the
    %   descending branch and handed back to restart a reversal curve at
    %   that point, so a reversal curve costs only its own points rather
    %   than a fresh descent from saturation.
    %
    %   Example:
    %       sim = FORCSimulation(1.0,-1.0,"Stride",5);
    %       sim.run();
    %       dist = sim.distribution(3);   % struct with Hr, H, rho, Hc, Hu
    %
    %   See also simulateHysteresisCurve.

    properties (SetAccess = private)
        Hmax (1,1) double
        Hmin (1,1) double
        SolveFcn

        HMajor (1,:) double % fields of the major descending branch
        MMajor (1,:) double % magnetization along the major descending branch

        % filled in by run(): the family of reversal curves on a common
        % grid, M(i,j) measured at field H(j) on the curve reversed at
        % Hr(i). Entries with H(j) < Hr(i) are NaN -- a reversal curve
        % says nothing about fields below the one it was reversed at.
        Hr (1,:) double
        H (1,:) double
        M double
    end

    properties (Access = private)
        history (1,:) cell      % one opaque solver state per HMajor step, plus the initial one
        reversalIndices (1,:) double
        ascendingIndices (1,:) cell
        Verbose (1,1) logical
    end

    methods
        function obj = FORCSimulation(Hmax,Hmin,options)
            %FORCSIMULATION Set up a reversal-curve study.
            %
            %   SIM = FORCSIMULATION(HMAX,HMIN) sweeps the field from
            %   HMAX down to HMIN to build the major descending branch,
            %   then chooses which points of it become reversal fields.
            %   Nothing is measured yet -- call SIM.RUN() for that.
            %
            %   Options:
            %       NumPoints         - number of field steps on the
            %                          descending branch (default 201)
            %       Stride            - use every STRIDE-th point of the
            %                          descending branch as a reversal
            %                          field (default 5); 1 uses them all
            %       MeasurementStride - measure on a finer grid than the
            %                          reversal fields are spaced on;
            %                          must divide Stride. Default keeps
            %                          the lattice square, which is what
            %                          the smoothing factor of the
            %                          distribution assumes.
            %       HrMin, HrMax      - restrict the reversal fields to
            %                          this window (default: the whole
            %                          sweep)
            %       SolveFcn          - the single-point magnetization
            %                          solver, called as
            %                          [M,state] = SolveFcn(H,Hmax,Hmin,state)
            %                          (default @simulateHysteresisCurve)
            %       Verbose           - report progress while running
            %                          (default true)
            arguments (Input)
                Hmax (1,1) double {mustBeFinite, mustBeReal}
                Hmin (1,1) double {mustBeFinite, mustBeReal}
                options.NumPoints (1,1) double {mustBePositive, mustBeInteger} = 201
                options.Stride (1,1) double {mustBePositive, mustBeInteger} = 5
                options.MeasurementStride double {mustBePositive, mustBeInteger} = []
                options.HrMin (1,1) double {mustBeReal} = -Inf
                options.HrMax (1,1) double {mustBeReal} = Inf
                options.SolveFcn = @simulateHysteresisCurve
                options.Verbose (1,1) logical = true
            end

            if Hmin >= Hmax
                error("FORCSimulation:BadWindow", ...
                      "Hmin (%g) must be below Hmax (%g).",Hmin,Hmax);
            end

            obj.Hmax = Hmax;
            obj.Hmin = Hmin;
            obj.SolveFcn = options.SolveFcn;
            obj.Verbose = options.Verbose;

            obj.simulateMajorBranch(options.NumPoints);
            obj.chooseReversalFields(options.Stride,options.MeasurementStride, ...
                options.HrMin,options.HrMax);
        end

        function run(obj)
            %RUN Measure every reversal curve and assemble the diagram.
            %
            %   SIM.RUN() fills in SIM.Hr, SIM.H and SIM.M. A curve whose
            %   field steps raise an error is reported and skipped
            %   instead of aborting the whole study.
            nBranches = numel(obj.reversalIndices);

            % the deepest reversal curve -- reversed latest, at the
            % smallest Hr -- carries the whole measurement grid; every
            % other curve measures the tail of it, since reversing later
            % means starting further up the same field ladder
            longest = obj.ascendingIndices{end};
            gridH = [obj.HMajor(obj.reversalIndices(end)) obj.HMajor(longest)];
            nH = numel(gridH);

            obj.Hr = obj.HMajor(obj.reversalIndices);
            obj.H = gridH;
            obj.M = NaN(nBranches,nH);

            nFailed = 0;
            for i = 1:nBranches
                k = obj.reversalIndices(i);
                measured = obj.ascendingIndices{i};
                obj.report("  reversal %d of %d (H = %.4g): %d points\n", ...
                    i,nBranches,obj.HMajor(k),numel(measured));

                try
                    ascendingM = obj.measureBranch(k,measured);
                catch ME
                    nFailed = nFailed + 1;
                    fprintf(2,"  reversal at H = %.4g failed: %s\n", ...
                        obj.HMajor(k),ME.message);
                    continue
                end

                % the value at the reversal field is not measured: it is
                % the one the major descending branch already has there
                branchM = [obj.MMajor(k) reshape(ascendingM,1,[])];
                nPoints = numel(branchM);
                obj.M(i,(nH - nPoints + 1):nH) = branchM;
            end

            obj.report("Study finished: %d reversal curves, %d failed.\n", ...
                nBranches - nFailed,nFailed);

            if nFailed == nBranches
                error("FORCSimulation:NoCurves", ...
                      "Every reversal curve failed; there is no diagram to build.");
            end
        end

        function d = distribution(obj,smoothingFactor)
            %DISTRIBUTION The FORC distribution of the measured family.
            %
            %   D = SIM.DISTRIBUTION(SF) fits
            %
            %       M ~ a1 + a2*u + a3*u^2 + a4*v + a5*v^2 + a6*u*v
            %
            %   by least squares over the (2*SF+1)-by-(2*SF+1) block of
            %   grid points around each point of the diagram, with u and
            %   v the offsets in Hr and H, and reads the mixed second
            %   derivative rho = -1/2 * a6 off the fit (Pike's method).
            %   Noisy families cannot be differentiated twice by finite
            %   differences alone, which is why this is a fit rather than
            %   a difference of neighbouring points.
            %
            %   SF is both the smoothing width and the resolution limit:
            %   a larger SF gives a cleaner but blunter distribution. rho
            %   is NaN wherever the window does not fit entirely inside
            %   the measured wedge H >= Hr (the diagonal and outer edges).
            %
            %   D is a plain struct with fields Hr, H, rho, Hc, Hu and
            %   SmoothingFactor. Hc = (H-Hr)/2 and Hu = (H+Hr)/2 are
            %   returned as matrices, one entry per (Hr,H) grid point,
            %   since that rotated grid is curvilinear.
            arguments (Input)
                obj (1,1) FORCSimulation
                smoothingFactor (1,1) double {mustBePositive, mustBeInteger} = 2
            end

            if isempty(obj.M)
                error("FORCSimulation:NotRun", ...
                      "Call run() before asking for the distribution.");
            end

            sf = smoothingFactor;
            window = 2*sf + 1;
            nHr = numel(obj.Hr);
            nH = numel(obj.H);
            if nHr < window || nH < window
                error("FORCSimulation:GridTooSmall", ...
                      "A smoothing factor of %d needs at least %d points along each axis; the grid is %d by %d.", ...
                      sf,window,nHr,nH);
            end

            rho = NaN(nHr,nH);
            for i = (1+sf):(nHr-sf)
                rowsIdx = (i-sf):(i+sf);
                for j = (1+sf):(nH-sf)
                    colsIdx = (j-sf):(j+sf);

                    block = obj.M(rowsIdx,colsIdx);
                    if any(isnan(block(:)))
                        % the window reaches outside the measured wedge
                        continue
                    end

                    [V,U] = meshgrid(obj.H(colsIdx) - obj.H(j),obj.Hr(rowsIdx) - obj.Hr(i));
                    u = U(:);
                    v = V(:);

                    su = max(abs(u));
                    sv = max(abs(v));
                    if su == 0 || sv == 0
                        continue
                    end
                    u = u/su;
                    v = v/sv;

                    A = [ones(numel(u),1), u, u.^2, v, v.^2, u.*v];
                    coefficients = A\block(:);

                    rho(i,j) = -0.5*coefficients(6)/(su*sv);
                end
            end

            [HH,HHr] = meshgrid(obj.H,obj.Hr);
            d = struct("Hr",obj.Hr,"H",obj.H,"rho",rho, ...
                "Hc",(HH - HHr)/2,"Hu",(HH + HHr)/2, ...
                "SmoothingFactor",sf);
        end
    end

    methods (Access = private)
        function simulateMajorBranch(obj,numPoints)
            %SIMULATEMAJORBRANCH Sweep from Hmax down to Hmin, recording state.
            %   One call to SolveFcn per field step; the solver state
            %   after every step is kept so a reversal curve can restart
            %   from any point of this branch instead of re-sweeping from
            %   saturation.
            obj.HMajor = linspace(obj.Hmax,obj.Hmin,numPoints);
            obj.MMajor = zeros(1,numPoints);
            obj.history = cell(1,numPoints + 1);

            state = [];
            obj.history{1} = state;
            for i = 1:numPoints
                [obj.MMajor(i),state] = obj.SolveFcn(obj.HMajor(i),obj.Hmax,obj.Hmin,state);
                obj.history{i+1} = state;
            end

            obj.report("Descending branch: %d field points from %.4g to %.4g.\n", ...
                numPoints,obj.Hmax,obj.Hmin);
        end

        function chooseReversalFields(obj,stride,measurementStride,hrMin,hrMax)
            %CHOOSEREVERSALFIELDS Pick which points of HMajor start a reversal curve.
            %
            %   A reversal curve starts at Hr = HMajor(k) and climbs back
            %   to HMajor(1), measured every MEASUREMENTSTRIDE-th master
            %   point on the way, so every curve lands on one common
            %   field grid. The first master point is saturation, where
            %   reversing would give a curve with no measurement points,
            %   so candidates start one stride below it.
            n = numel(obj.HMajor);
            if isempty(measurementStride)
                measurementStride = stride;
            end
            if mod(stride,measurementStride) ~= 0
                error("FORCSimulation:IncommensurateStrides", ...
                      "A measurement stride of %d does not divide the reversal stride of %d, so the curves would not share one field grid.", ...
                      measurementStride,stride);
            end
            if hrMin > hrMax
                error("FORCSimulation:EmptyWindow", ...
                      "HrMin (%g) is above HrMax (%g), so no reversal field can be selected.", ...
                      hrMin,hrMax);
            end

            candidates = (1 + stride):stride:n;
            inWindow = obj.HMajor(candidates) >= hrMin & obj.HMajor(candidates) <= hrMax;
            indices = candidates(inWindow);

            if isempty(indices)
                error("FORCSimulation:EmptyPlan", ...
                      "No reversal field left after applying stride %d and the window [%g, %g] to %d master points.", ...
                      stride,hrMin,hrMax,n);
            end

            ascending = cell(1,numel(indices));
            for i = 1:numel(indices)
                ascending{i} = (indices(i) - measurementStride):-measurementStride:1;
            end

            obj.reversalIndices = indices;
            obj.ascendingIndices = ascending;

            nSolves = sum(cellfun(@numel,ascending));
            obj.report("%d reversal curves, %d equilibrium solves (the descending branch was %d).\n", ...
                numel(indices),nSolves,n);
        end

        function ascendingM = measureBranch(obj,k,measured)
            %MEASUREBRANCH Restore the state at HMajor(k) and sweep back up.
            %   History row 1 is the state before any field step, so the
            %   state after the k-th field is history{k+1}.
            state = obj.history{k+1};
            Hascending = obj.HMajor(measured);
            ascendingM = zeros(1,numel(Hascending));
            for j = 1:numel(Hascending)
                [ascendingM(j),state] = obj.SolveFcn(Hascending(j),obj.Hmax,obj.Hmin,state);
            end
        end

        function report(obj,format,varargin)
            %REPORT Print progress unless the study was asked to stay quiet.
            if obj.Verbose
                fprintf(format,varargin{:});
            end
        end
    end
end
