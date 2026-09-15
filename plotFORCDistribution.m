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
    options.OutputPrefix (1,1) string = ""
    options.MarkerHcHu (:,2) double = zeros(0,2)
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

d = forcDistribution(Hr,H,M,options.SmoothingFactor);
peak = max(abs(d.rho),[],"all","omitnan");

[f,ax] = createFigure(FONTSIZE);
f.Position(3:4) = FIGURE_SIZE;

% The (Hc,Hu) grid is the (Hr,H) lattice rotated by 45 degrees, so it is
% curvilinear and is passed as matrices. Interpolating it onto a straight
% grid first would blur the ridge this figure is drawn to show.
contourf(ax,d.Hc,d.Hu,d.rho,CONTOUR_LEVELS,"LineColor","none","HandleVisibility","off")
colormap(ax,divergingColormap)
clim(ax,[-peak,peak])

yline(ax,0,"LineStyle",":","Color","k","LineWidth",0.8,"HandleVisibility","off")

if ~isempty(options.MarkerHcHu)
    plot(ax,options.MarkerHcHu(:,1),options.MarkerHcHu(:,2), ...
        "Marker","x","MarkerSize",14,"LineWidth",2.5,"LineStyle","none", ...
        "Color","k","DisplayName","Theoretical")
    legend(ax,"Location","northeast","FontSize",0.6*FONTSIZE)
end

xlabel(ax,"Coercivity coordinate, H_c")
ylabel(ax,"Interaction coordinate, H_u")

cb = colorbar(ax);
cb.Label.String = "FORC distribution, \rho";

title(ax,sprintf("Smoothing factor SF = %d",options.SmoothingFactor), ...
    "FontSize",0.75*FONTSIZE)
hold(ax,"off")

savePlot(f,options.OutputPrefix + "figFORCDistribution.pdf")

end
