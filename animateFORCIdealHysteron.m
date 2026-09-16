function animateFORCIdealHysteron(outputFile,options)
%ANIMATEFORCIDEALHYSTERON Animate ideal hysteron switches becoming FORC points.
%
%   ANIMATEFORCIDEALHYSTERON(OUTPUTFILE) reads the JSON written by
%   runFORCDemoIdealHysteron.m (one particle) or
%   runFORCDemoTwoHysterons.m (several independent particles) and
%   renders a two-panel animated GIF:
%
%       left  -- the applied field sweeping down the major branch, then
%                 each reversal curve traced back up to saturation
%       right -- the FORC distribution in (Hc,Hu) coordinates
%
%   The right panel stays blank while curves are being traced: Pike's
%   fit needs a full (2*SF+1)-by-(2*SF+1) neighbourhood of MEASURED
%   points, so there is no honest way to show a partial distribution
%   part-way through the family -- only once every curve selected for
%   display has been traced does the fit run and the point/ridge
%   appear, marked against the theoretical (Hc,Hu) this demo predicts.
%
%   The switching field Hsw is marked with a thick dashed line in both
%   panels throughout: at H = +-Hsw on the left (where the particle
%   actually switches), and at Hc = Hsw on the right (where that
%   switch shows up in the distribution) -- so it is visually obvious
%   that the two are the same number, not just similar-looking peaks.
%
%   OUTPUTFILE defaults to "forc_demo_ideal_hysteron_output.json" next
%   to this script. The distribution itself is always computed from
%   the FULL measured family in that file -- MaxCurves only controls
%   how many reversal curves are drawn in the animation, not what the
%   fit is based on.
%
%   Options:
%       OutputFile      - GIF filename, written into plots/
%                         (default "idealHysteronAnimation.gif")
%       FrameRate       - playback frame rate (frames per second); lower
%                         plays back slower without changing the
%                         animation itself (default 15)
%       MaxCurves       - number of reversal curves to animate, evenly
%                         sampled like plotFORCFamily's MaxCurves
%                         (default Inf, i.e. every simulated reversal
%                         curve; pass a finite value to thin the
%                         animation and shorten the video)
%       SmoothingFactor - passed to forcDistribution for the final
%                         reveal (default 2)
%       DescentFrames   - frames used to trace each curve's descent to
%                         its reversal field (default 8)
%       AscentFrames    - frames used to trace each curve's ascent back
%                         to saturation (default 12)
%       HoldFramesAtEnd - frames the finished distribution is held for,
%                         so it survives being paused on a slide
%                         (default 40)
%
%   See also runFORCDemoIdealHysteron, runFORCDemoTwoHysterons,
%   plotFORCFamily, plotFORCDistribution, forcDistribution.

arguments (Input)
    outputFile (1,1) string = ""
    options.OutputFile (1,1) string = "idealHysteronAnimation.gif"
    options.FrameRate (1,1) double {mustBePositive} = 15
    options.MaxCurves (1,1) double {mustBePositive} = Inf
    options.SmoothingFactor (1,1) double {mustBePositive, mustBeInteger} = 2
    options.DescentFrames (1,1) double {mustBePositive, mustBeInteger} = 8
    options.AscentFrames (1,1) double {mustBePositive, mustBeInteger} = 12
    options.HoldFramesAtEnd (1,1) double {mustBeNonnegative, mustBeInteger} = 40
end

if outputFile == ""
    outputFile = fullfile(fileparts(mfilename("fullpath")),"forc_demo_ideal_hysteron_output.json");
end

FONTSIZE = 20;
MAJOR_COLOR = [0.75 0.75 0.75];
RECORDED_COLOR = [0.35 0.35 0.35];
HIGHLIGHT_COLOR = [0.85 0.10 0.10];
SWITCHING_COLOR = [0.90 0.60 0.00];
CONTOUR_LEVELS = 24;

raw = jsondecode(fileread(outputFile));
Hr = raw.diagram.Hr(:)';
H = raw.diagram.H(:)';
M = raw.diagram.M;
Hmajor = raw.majorBranch.H(:)';
Mmajor = raw.majorBranch.M(:)';

if isfield(raw,"theoretical")
    markerHcHu = [raw.theoretical.Hc, raw.theoretical.Hu];
else
    markerHcHu = zeros(0,2);
end

if isfield(raw,"parameters") && isfield(raw.parameters,"Hsw")
    Hsw = raw.parameters.Hsw(:)';
elseif ~isempty(markerHcHu)
    Hsw = markerHcHu(:,1)';
else
    Hsw = [];
end

% Computed once up front: MaxCurves only thins which curves are drawn,
% never what the fit sees, so the distribution the final reveal shows
% is exactly the one plotFORCDistribution would draw from this file.
d = forcDistribution(Hr,H,M,options.SmoothingFactor);
peak = max(abs(d.rho),[],"all","omitnan");

nCurves = min(options.MaxCurves,numel(Hr));
rows = round(linspace(1,numel(Hr),nCurves));

plotsDir = fullfile(fileparts(mfilename("fullpath")),"plots");
if ~isfolder(plotsDir)
    mkdir(plotsDir)
end
[~,baseName] = fileparts(options.OutputFile);
gifPath = fullfile(plotsDir,baseName + ".gif");

f = figure("Renderer","painters");
fontsize(f,FONTSIZE,"points")
f.Position(3:4) = [1500 700];
tl = tiledlayout(f,1,2,"TileSpacing","compact","Padding","compact");

axLeft = nexttile(tl);
hold(axLeft,"on"); grid(axLeft,"on"); axis(axLeft,"square")
xlim(axLeft,[min(Hmajor) max(Hmajor)])
ylim(axLeft,[-1.1 1.1])
xlabel(axLeft,"Applied field, H")
ylabel(axLeft,"Magnetization, M")
title(axLeft,"Field sweep and reversal curves","FontSize",0.8*FONTSIZE)
plot(axLeft,Hmajor,Mmajor,"Color",MAJOR_COLOR,"LineWidth",2,"HandleVisibility","off")
if ~isempty(Hsw)
    for i = 1:numel(Hsw)
        if numel(Hsw) == 1
            label = "H_{sw} (switching field)";
        else
            label = sprintf("H_{sw,%d} = %.2g",i,Hsw(i));
        end
        xline(axLeft,Hsw(i),"LineStyle","--","Color",SWITCHING_COLOR,"LineWidth",3, ...
            "DisplayName",label)
        xline(axLeft,-Hsw(i),"LineStyle","--","Color",SWITCHING_COLOR,"LineWidth",3, ...
            "HandleVisibility","off")
    end
    legend(axLeft,"Location","east","FontSize",0.6*FONTSIZE)
end
posMarker = plot(axLeft,NaN,NaN,"o","MarkerSize",10, ...
    "MarkerFaceColor",HIGHLIGHT_COLOR,"MarkerEdgeColor","k","HandleVisibility","off");
traceLine = animatedline(axLeft,"Color",HIGHLIGHT_COLOR,"LineWidth",2.2,"HandleVisibility","off");

axRight = nexttile(tl);
hold(axRight,"on"); grid(axRight,"on"); axis(axRight,"square")
xlim(axRight,[min(d.Hc,[],"all") max(d.Hc,[],"all")])
ylim(axRight,[min(d.Hu,[],"all") max(d.Hu,[],"all")])
xlabel(axRight,"Coercivity coordinate, H_c")
ylabel(axRight,"Interaction coordinate, H_u")
title(axRight,"FORC distribution (revealed once the family is complete)", ...
    "FontSize",0.7*FONTSIZE)
if ~isempty(Hsw)
    for i = 1:numel(Hsw)
        if numel(Hsw) == 1
            label = "H_c = H_{sw}";
        else
            label = sprintf("H_c = H_{sw,%d} = %.2g",i,Hsw(i));
        end
        xline(axRight,Hsw(i),"LineStyle","--","Color",SWITCHING_COLOR,"LineWidth",3, ...
            "DisplayName",label)
    end
end
caption = text(axRight,mean(xlim(axRight)),mean(ylim(axRight)), ...
    sprintf("curve 0 of %d",nCurves), ...
    "HorizontalAlignment","center","FontSize",FONTSIZE);

gifFrameIndex = 0;
frameSize = [];

for k = 1:nCurves
    row = rows(k);

    if k == 1
        startH = max(Hmajor);
    else
        startH = Hr(rows(k-1));
    end
    [~,idxStart] = min(abs(Hmajor - startH));
    [~,idxEnd] = min(abs(Hmajor - Hr(row)));

    descH = Hmajor(idxStart:idxEnd);
    descM = Mmajor(idxStart:idxEnd);
    descSel = subsampleIndices_(numel(descH),options.DescentFrames);
    for idx = descSel
        set(posMarker,"XData",descH(idx),"YData",descM(idx));
        addpoints(traceLine,descH(idx),descM(idx));
        [gifFrameIndex,frameSize] = emitFrame_(f,gifPath,gifFrameIndex,frameSize,options.FrameRate);
    end

    ascCols = find(~isnan(M(row,:)));
    ascH = H(ascCols);
    ascM = M(row,ascCols);
    ascSel = subsampleIndices_(numel(ascH),options.AscentFrames);
    for idx = ascSel
        set(posMarker,"XData",ascH(idx),"YData",ascM(idx));
        addpoints(traceLine,ascH(idx),ascM(idx));
        [gifFrameIndex,frameSize] = emitFrame_(f,gifPath,gifFrameIndex,frameSize,options.FrameRate);
    end

    plot(axLeft,ascH,ascM,"Color",RECORDED_COLOR,"LineWidth",1.1,"HandleVisibility","off")
    clearpoints(traceLine)

    set(caption,"String",sprintf("curve %d of %d",k,nCurves))
    [gifFrameIndex,frameSize] = emitFrame_(f,gifPath,gifFrameIndex,frameSize,options.FrameRate);
end

delete(caption)
contourf(axRight,d.Hc,d.Hu,d.rho,CONTOUR_LEVELS,"LineColor","none","HandleVisibility","off")
colormap(axRight,divergingColormap)
clim(axRight,[-peak,peak])
yline(axRight,0,"LineStyle",":","Color","k","LineWidth",0.8,"HandleVisibility","off")
if ~isempty(markerHcHu)
    plot(axRight,markerHcHu(:,1),markerHcHu(:,2), ...
        "Marker","x","MarkerSize",14,"LineWidth",2.5,"LineStyle","none","Color","k", ...
        "DisplayName","Theoretical (H_{sw},0)")
end
if ~isempty(Hsw) || ~isempty(markerHcHu)
    legend(axRight,"Location","northeast","FontSize",0.6*FONTSIZE)
end
cb = colorbar(axRight);
cb.Label.String = "FORC distribution, \rho";
title(axRight,sprintf("FORC distribution, SF = %d",options.SmoothingFactor), ...
    "FontSize",0.7*FONTSIZE)

for i = 1:options.HoldFramesAtEnd
    [gifFrameIndex,frameSize] = emitFrame_(f,gifPath,gifFrameIndex,frameSize,options.FrameRate);
end

hold(axLeft,"off")
hold(axRight,"off")

savePlot(f,baseName + "_finalFrame.pdf")

fprintf("Wrote %s\n",gifPath);

end


function idx = subsampleIndices_(n,maxFrames)
%SUBSAMPLEINDICES_ Evenly spaced indices 1:n, capped at maxFrames, always including both ends.
idx = unique(round(linspace(1,n,min(maxFrames,n))));
end


function [frameIndex,frameSize] = emitFrame_(f,gifPath,frameIndex,frameSize,frameRate)
%EMITFRAME_ Capture the current figure and append it to the GIF.
%
%   GETFRAME's pixel size can still drift if the figure's content
%   bounding box changes (e.g. the colorbar/legend added partway
%   through this animation), and every GIF frame must match the first
%   one it was given. FRAMESIZE is fixed from this function's first call
%   and every later frame is resized to match.
drawnow
frame = getframe(f);
img = frame.cdata;
if isempty(frameSize)
    frameSize = [size(img,1) size(img,2)];
elseif ~isequal([size(img,1) size(img,2)],frameSize)
    img = imresize(img,frameSize);
end
[indexed,map] = rgb2ind(img,256);
if frameIndex == 0
    imwrite(indexed,map,gifPath,"gif","LoopCount",Inf,"DelayTime",1/max(frameRate,1));
else
    imwrite(indexed,map,gifPath,"gif","WriteMode","append","DelayTime",1/max(frameRate,1));
end
frameIndex = frameIndex + 1;
end
