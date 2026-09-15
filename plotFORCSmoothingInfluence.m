function plotFORCSmoothingInfluence(outputFile,options)
%PLOTFORCSMOOTHINGINFLUENCE What the smoothing factor decides, measured.
%
%   PLOTFORCSMOOTHINGINFLUENCE(OUTPUTFILE) draws four panels over one
%   run's diagram, saved as figFORCSmoothingInfluenceA.pdf through
%   figFORCSmoothingInfluenceD.pdf:
%
%       A  the lowest Hc at which rho is defined, against SF*dH
%       B  the peak rho
%       C  the ridge centre <Hc>
%       D  the ridge width in Hu
%
%   Panels C and D each carry two series, and the difference between them
%   is the point of the figure. Raising SF deletes the lowest-Hc columns
%   of the diagram, so a moment summed over the native support is summed
%   over a different region at every SF. Repeating the sum on the support
%   the largest SF leaves -- which is contained in all the others --
%   removes that effect and leaves only the smoothing.
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
MARKER_SIZE = 10;
LINE_WIDTH = 1.5;
NATIVE_COLOR = [0.85 0.33 0.10];
FIXED_COLOR = [0.15 0.55 0.20];
PREDICTION_COLOR = [0.45 0.45 0.45];

raw = jsondecode(fileread(outputFile));
Hr = raw.diagram.Hr(:)';
H = raw.diagram.H(:)';
M = raw.diagram.M;
dH = abs(median(diff(H)));

sfValues = options.SFValues;
rho = cell(1,numel(sfValues));
for isf = 1:numel(sfValues)
    rho{isf} = forcDistribution(Hr,H,M,sfValues(isf));
end

% The support shrinks as SF grows and the largest one is contained in
% all the others, so it is the only region on which distributions from
% different SF can be compared without the comparison being an artefact
% of the region summed over.
commonSupport = isfinite(rho{end}.rho);

lowestHc = arrayfun(@(isf) lowestHc_(rho{isf}),1:numel(sfValues));
peakRho = arrayfun(@(isf) max(rho{isf}.rho,[],"all","omitnan"),1:numel(sfValues));
[nativeCentre,nativeWidth] = deal(zeros(1,numel(sfValues)));
[fixedCentre,fixedWidth] = deal(zeros(1,numel(sfValues)));
for isf = 1:numel(sfValues)
    [nativeCentre(isf),nativeWidth(isf)] = ridgeMoments(rho{isf},[]);
    [fixedCentre(isf),fixedWidth(isf)] = ridgeMoments(rho{isf},commonSupport);
end

% Panel A: the cut is the smoothing window meeting the H = Hr boundary,
% so it should sit exactly on SF*dH -- and a cut that tracks SF is not an
% absence of low-coercivity switching.
[fA,axA] = createFigure(FONTSIZE);
q(1) = plot(axA,sfValues,sfValues*dH, ...
    "Color",PREDICTION_COLOR,"LineWidth",LINE_WIDTH,"DisplayName","SF x dH");
q(2) = plot(axA,sfValues,lowestHc, ...
    "LineStyle","none","Marker","o","MarkerSize",MARKER_SIZE, ...
    "MarkerFaceColor","b","MarkerEdgeColor","k","DisplayName","Measured");
ylabel(axA,"Lowest defined H_c")
title(axA,"The cut is exactly SF x dH","FontSize",0.8*FONTSIZE)
legend(axA,q,"Location","northwest","FontSize",0.8*FONTSIZE)

% Panel B: peak rho is a single-cell statistic on a noisy second
% derivative, and on top of that it depends on SF by an order of
% magnitude. It should not be quoted at all.
[fB,axB] = createFigure(FONTSIZE);
plot(axB,sfValues,peakRho, ...
    "Color","r","LineWidth",LINE_WIDTH,"Marker","o", ...
    "MarkerSize",MARKER_SIZE,"MarkerFaceColor","r")
set(axB,"YScale","log")
ylabel(axB,"Peak \rho")
title(axB,sprintf("Amplitude falls %.0fx over the sweep", ...
    peakRho(1)/peakRho(end)),"FontSize",0.8*FONTSIZE)

[fC,axC] = createFigure(FONTSIZE);
plotSupportPair(axC,sfValues,nativeCentre,fixedCentre,sfValues(end), ...
    NATIVE_COLOR,FIXED_COLOR,MARKER_SIZE,LINE_WIDTH,FONTSIZE)
ylabel(axC,"Ridge centre, <H_c>")
title(axC,"The coercivity drift is domain truncation","FontSize",0.8*FONTSIZE)

[fD,axD] = createFigure(FONTSIZE);
plotSupportPair(axD,sfValues,nativeWidth,fixedWidth,sfValues(end), ...
    NATIVE_COLOR,FIXED_COLOR,MARKER_SIZE,LINE_WIDTH,FONTSIZE)
ylabel(axD,"Ridge width in H_u")
title(axD,"The interaction width drifts even at fixed domain","FontSize",0.8*FONTSIZE)

for ax = [axA,axB,axC,axD]
    xlabel(ax,"Smoothing factor, SF")
    xticks(ax,sfValues)
    xlim(ax,[sfValues(1)-0.5,sfValues(end)+0.5])
    hold(ax,"off")
end

savePlot(fA,options.OutputPrefix + "figFORCSmoothingInfluenceA.pdf")
savePlot(fB,options.OutputPrefix + "figFORCSmoothingInfluenceB.pdf")
savePlot(fC,options.OutputPrefix + "figFORCSmoothingInfluenceC.pdf")
savePlot(fD,options.OutputPrefix + "figFORCSmoothingInfluenceD.pdf")

end


function value = lowestHc_(d)
Hc = d.Hc;
value = min(Hc(isfinite(d.rho)));
end


function [centreHc,widthHu] = ridgeMoments(d,mask)
% Weighted centre in Hc and full width at half maximum in Hu of the ridge.
%
% Only the positive part is summed. rho weights the moments, so including
% the negative lobe would subtract from the very density the moments are
% meant to describe rather than broaden it.
keep = isfinite(d.rho) & d.rho > 0;
if ~isempty(mask)
    keep = keep & mask;
end
weights = d.rho(keep);
Hc = d.Hc;
Hu = d.Hu;
centreHc = sum(weights.*Hc(keep))/sum(weights);
centreHu = sum(weights.*Hu(keep))/sum(weights);
variance = sum(weights.*(Hu(keep) - centreHu).^2)/sum(weights);
widthHu = 2*sqrt(2*log(2))*sqrt(variance);
end


function plotSupportPair(ax,sfValues,nativeValue,fixedValue,largestSF, ...
    nativeColor,fixedColor,markerSize,lineWidth,fontSize)
% One moment measured both ways, so the domain effect can be read off.
p(1) = plot(ax,sfValues,nativeValue, ...
    "Color",nativeColor,"LineWidth",lineWidth,"Marker","o", ...
    "MarkerSize",markerSize,"MarkerFaceColor",nativeColor, ...
    "DisplayName","Native support");
p(2) = plot(ax,sfValues,fixedValue, ...
    "Color",fixedColor,"LineStyle","--","LineWidth",lineWidth,"Marker","s", ...
    "MarkerSize",markerSize,"MarkerFaceColor",fixedColor, ...
    "DisplayName",sprintf("Support of SF = %d",largestSF));
legend(ax,p,"Location","best","FontSize",0.7*fontSize);
end
