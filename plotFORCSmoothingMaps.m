function plotFORCSmoothingMaps(outputFile,options)
%PLOTFORCSMOOTHINGMAPS One FORC distribution at several smoothing factors.
%
%   PLOTFORCSMOOTHINGMAPS(OUTPUTFILE) draws rho at each smoothing factor
%   in SFValues, one figure per factor, saved as
%   figFORCSmoothingMapSF<k>.pdf.
%
%   The smoothing factor is usually presented as a cosmetic choice. It is
%   not. rho at a grid point is a least squares fit over the surrounding
%   (2*SF+1)^2 block, and that block must lie entirely inside the
%   measured wedge H >= Hr, so raising SF does three things at once: it
%   smooths the texture, it collapses the amplitude, and it deletes the
%   lowest-Hc columns of the diagram altogether. The dashed line marks
%   that last effect, at Hc = SF*dH.
%
%   Each panel is normalised to its own peak, and the peak is stated in
%   the title instead: the amplitude falls sharply across the sweep, and
%   a shared colour scale would leave the high-SF panels blank.
%
%   OUTPUTFILE defaults to "forc_demo_output.json" next to this script.

arguments (Input)
    outputFile (1,1) string = ""
    options.SFValues (1,:) double {mustBePositive, mustBeInteger} = 1:5
    options.OutputPrefix (1,1) string = ""
end

if outputFile == ""
    outputFile = fullfile(fileparts(mfilename("fullpath")),"forc_demo_output.json");
end

FONTSIZE = 24;
CONTOUR_LEVELS = 24;

raw = jsondecode(fileread(outputFile));
Hr = raw.diagram.Hr(:)';
H = raw.diagram.H(:)';
M = raw.diagram.M;

% The lattice is square by construction, so one spacing describes both
% axes and the cut can be quoted as a single field.
dH = abs(median(diff(H)));

for sf = options.SFValues
    d = forcDistribution(Hr,H,M,sf);
    peak = max(abs(d.rho),[],"all","omitnan");

    [f,ax] = createFigure(FONTSIZE);

    contourf(ax,d.Hc,d.Hu,d.rho/peak,CONTOUR_LEVELS,"LineColor","none")
    colormap(ax,divergingColormap)
    clim(ax,[-1,1])

    % where the fitting window first fits inside the measured wedge
    xline(ax,sf*dH,"LineStyle","--","Color","k","LineWidth",1.4)
    yline(ax,0,"LineStyle",":","Color","k","LineWidth",0.8)
    xlabel(ax,"Coercivity coordinate, H_c")
    ylabel(ax,"Interaction coordinate, H_u")

    cb = colorbar(ax);
    cb.Label.String = "\rho / \rho_{peak}";

    title(ax,sprintf("SF = %d (%d x %d), cut %.3g, peak %.3g", ...
        sf,2*sf+1,2*sf+1,sf*dH,peak),"FontSize",0.7*FONTSIZE)
    hold(ax,"off")

    savePlot(f,options.OutputPrefix + sprintf("figFORCSmoothingMapSF%d.pdf",sf))
end

end
