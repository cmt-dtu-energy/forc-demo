function plotFORCFamily(outputFile,options)
%PLOTFORCFAMILY The family of first-order reversal curves of one run.
%
%   PLOTFORCFAMILY(OUTPUTFILE) plots every reversal curve measured by
%   runFORCDemo, with the major descending branch they all leave off
%   drawn on top. Nothing here is a second simulation: each curve
%   branches off the descending branch at its own reversal field, so the
%   family is the raw measurement the FORC distribution is computed
%   from, and this figure is what shows that the two are the same data.
%
%   OUTPUTFILE defaults to "forc_demo_output.json" next to this script,
%   i.e. the default output of runFORCDemo.

arguments (Input)
    outputFile (1,1) string = ""
    options.MaxCurves (1,1) double {mustBePositive, mustBeInteger} = 10
end

if outputFile == ""
    outputFile = fullfile(fileparts(mfilename("fullpath")),"forc_demo_output.json");
end

FONTSIZE = 22;
CURVE_COLOR = [0.35 0.35 0.35];
CURVE_LINE_WIDTH = 1;
CURVE_OPACITY = 0.7;
MAJOR_LINE_WIDTH = 2.0;

raw = jsondecode(fileread(outputFile));
Hr = raw.diagram.Hr(:)';
H = raw.diagram.H(:)';
M = raw.diagram.M;              % jsondecode reads the unmeasured wedge back as NaN
Hmajor = raw.majorBranch.H(:)';
Mmajor = raw.majorBranch.M(:)';

[f,ax] = createFigure(FONTSIZE);

% The whole family in one call, one line per row: the NaN wedge is
% dropped by plot itself, so every curve starts at its own reversal
% field without any of them being trimmed by hand.
nCurves = min(options.MaxCurves,numel(Hr));
rows = round(linspace(1,numel(Hr),nCurves));
h = plot(ax,H,M(rows,:)', ...
    "Color",CURVE_COLOR, ...
    "LineWidth",CURVE_LINE_WIDTH, ...
    "HandleVisibility","off");
for i = 1:numel(h)
    h(i).Color(4) = CURVE_OPACITY;    % fourth channel is the line opacity
end

% proxy so the whole family takes a single legend entry
p(1) = plot(ax,NaN,NaN, ...
    "Color",CURVE_COLOR, ...
    "LineWidth",CURVE_LINE_WIDTH, ...
    "DisplayName","Reversal curves");
p(2) = plot(ax,Hmajor,Mmajor, ...
    "Color","k", ...
    "LineWidth",MAJOR_LINE_WIDTH, ...
    "DisplayName","Major descending branch");

xline(ax,0,"LineStyle",":","Color","k","LineWidth",0.6, ...
    "HandleVisibility","off")
yline(ax,0,"LineStyle",":","Color","k","LineWidth",0.6, ...
    "HandleVisibility","off")
xlabel(ax,"Applied field, H")
ylabel(ax,"Magnetization, M")
legend(ax,p,"Location","northwest","FontSize",0.7*FONTSIZE);
title(ax,"FORC family (placeholder solver)","FontSize",0.75*FONTSIZE)
hold(ax,"off")

savePlot(f,"figFORCFamily.pdf")

end
