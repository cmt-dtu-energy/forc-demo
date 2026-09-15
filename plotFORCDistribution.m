function plotFORCDistribution(outputFile,options)
%PLOTFORCDISTRIBUTION The FORC distribution of one run.
%
%   PLOTFORCDISTRIBUTION(OUTPUTFILE) draws rho = -1/2 d2M/(dHr dH) in the
%   customary rotated coordinates
%
%       Hc = (H - Hr)/2      Hu = (H + Hr)/2
%
%   read as the joint density of coercivities and of the "interaction"
%   fields the sample's internal memory imposes: spread along Hc is the
%   coercivity spread, spread along Hu the interaction spread.
%
%   The colour scale is diverging and symmetric about zero: rho changes
%   sign, and a negative lobe is a different statement than a weak
%   positive one, so zero has to be the colour the eye reads as neutral.
%
%   OUTPUTFILE defaults to "forc_demo_output.json" next to this script.

arguments (Input)
    outputFile (1,1) string = ""
    options.SmoothingFactor (1,1) double {mustBePositive, mustBeInteger} = 3
end

if outputFile == ""
    outputFile = fullfile(fileparts(mfilename("fullpath")),"forc_demo_output.json");
end

FONTSIZE = 26;
CONTOUR_LEVELS = 24;
FIGURE_SIZE = [900 700];

raw = jsondecode(fileread(outputFile));
Hr = raw.diagram.Hr(:)';
H = raw.diagram.H(:)';
M = raw.diagram.M;

d = forcDistribution_(Hr,H,M,options.SmoothingFactor);
peak = max(abs(d.rho),[],"all","omitnan");

[f,ax] = createFigure(FONTSIZE);
f.Position(3:4) = FIGURE_SIZE;

% The (Hc,Hu) grid is the (Hr,H) lattice rotated by 45 degrees, so it is
% curvilinear and is passed as matrices. Interpolating it onto a straight
% grid first would blur the ridge this figure is drawn to show.
contourf(ax,d.Hc,d.Hu,d.rho,CONTOUR_LEVELS,"LineColor","none")
colormap(ax,divergingColormap)
clim(ax,[-peak,peak])

yline(ax,0,"LineStyle",":","Color","k","LineWidth",0.8)
xlabel(ax,"Coercivity coordinate, H_c")
ylabel(ax,"Interaction coordinate, H_u")

cb = colorbar(ax);
cb.Label.String = "FORC distribution, \rho";

title(ax,sprintf("Smoothing factor SF = %d",options.SmoothingFactor), ...
    "FontSize",0.75*FONTSIZE)
hold(ax,"off")

savePlot(f,"figFORCDistribution.pdf")

end


function d = forcDistribution_(Hr,H,M,smoothingFactor)
%FORCDISTRIBUTION_ Pike's mixed second derivative of a FORC family.
%
%   Fits M ~ a1 + a2*u + a3*u^2 + a4*v + a5*v^2 + a6*u*v by least squares
%   over the (2*SF+1)-by-(2*SF+1) neighbourhood of each grid point (u, v
%   the offsets in Hr, H) and reads the mixed derivative off as a6. See
%   FORCSimulation.distribution for the same fit, kept here as a small
%   local copy so this plotting script has no non-plotting dependency.
sf = smoothingFactor;
nHr = numel(Hr);
nH = numel(H);
rho = NaN(nHr,nH);
for i = (1+sf):(nHr-sf)
    rowsIdx = (i-sf):(i+sf);
    for j = (1+sf):(nH-sf)
        colsIdx = (j-sf):(j+sf);
        block = M(rowsIdx,colsIdx);
        if any(isnan(block(:)))
            continue
        end
        [V,U] = meshgrid(H(colsIdx) - H(j),Hr(rowsIdx) - Hr(i));
        u = U(:); v = V(:);
        su = max(abs(u)); sv = max(abs(v));
        if su == 0 || sv == 0
            continue
        end
        u = u/su; v = v/sv;
        A = [ones(numel(u),1), u, u.^2, v, v.^2, u.*v];
        coefficients = A\block(:);
        rho(i,j) = -0.5*coefficients(6)/(su*sv);
    end
end
[HH,HHr] = meshgrid(H,Hr);
d = struct("Hr",Hr,"H",H,"rho",rho,"Hc",(HH-HHr)/2,"Hu",(HH+HHr)/2, ...
    "SmoothingFactor",sf);
end
