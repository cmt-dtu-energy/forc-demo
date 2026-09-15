function animateFORCIdealHysteron(outputFile,options)
%ANIMATEFORCIDEALHYSTERON Animate ideal hysteron switches becoming FORC points.
%
%   ANIMATEFORCIDEALHYSTERON(OUTPUTFILE) reads the JSON written by
%   runFORCDemoIdealHysteron.m (one particle) or
%   runFORCDemoTwoHysterons.m (several independent particles) and
%   renders a two-panel video:
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
%       OutputFile      - video filename, written into plots/
%                         (default "idealHysteronAnimation.mp4")
%       FrameRate       - video frame rate; lower plays back slower
%                         without changing the animation itself
%                         (default 15)
%       MaxCurves       - number of reversal curves to animate, evenly
%                         sampled like plotFORCFamily's MaxCurves
%                         (default 12)
%       SmoothingFactor - passed to forcDistribution for the final
%                         reveal (default 2)
%       DescentFrames   - frames used to trace each curve's descent to
%                         its reversal field (default 8)
%       AscentFrames    - frames used to trace each curve's ascent back
%                         to saturation (default 12)
%       HoldFramesAtEnd - frames the finished distribution is held for,
%                         so it survives being paused on a slide
%                         (default 40)
%       AlsoWriteGif    - also write an animated GIF alongside the MP4,
%                         for README/chat use rather than slides
%                         (default false)
%
%   See also runFORCDemoIdealHysteron, runFORCDemoTwoHysterons,
%   plotFORCFamily, plotFORCDistribution, forcDistribution.

arguments (Input)
    outputFile (1,1) string = ""
    options.OutputFile (1,1) string = "idealHysteronAnimation.mp4"
    options.FrameRate (1,1) double {mustBePositive} = 15
    options.MaxCurves (1,1) double {mustBePositive, mustBeInteger} = 12
    options.SmoothingFactor (1,1) double {mustBePositive, mustBeInteger} = 2
    options.DescentFrames (1,1) double {mustBePositive, mustBeInteger} = 8
    options.AscentFrames (1,1) double {mustBePositive, mustBeInteger} = 12
    options.HoldFramesAtEnd (1,1) double {mustBeNonnegative, mustBeInteger} = 40
    options.AlsoWriteGif (1,1) logical = false
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
% Recorded as Motion JPEG in an AVI container, not straight to MP4: this
% repo's usual MPEG-4/H.264 VideoWriter profile was found to silently
% write a corrupted bitstream in this environment (confirmed with an
% independent decoder, ffmpeg, not just MATLAB's own VideoReader --
% verified frame-by-frame that the corruption is baked into the file,
% not an artifact of how it's read back). Motion JPEG held up cleanly
% across the same test. If `ffmpeg` is on the system path it's used
% afterwards to transcode this into the requested MP4 (smaller, and the
% format slide software actually expects); otherwise the AVI itself is
% the delivered output, so this never silently ships a broken video.
aviPath = fullfile(plotsDir,baseName + ".avi");
gifPath = fullfile(plotsDir,baseName + ".gif");

v = VideoWriter(aviPath,"Motion JPEG AVI");
v.FrameRate = options.FrameRate;
v.Quality = 90;
open(v);
closeVideo = onCleanup(@() closeIfOpen_(v));

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
        [v,gifFrameIndex,frameSize] = emitFrame_(v,f,gifPath,gifFrameIndex,frameSize,options.AlsoWriteGif);
    end

    ascCols = find(~isnan(M(row,:)));
    ascH = H(ascCols);
    ascM = M(row,ascCols);
    ascSel = subsampleIndices_(numel(ascH),options.AscentFrames);
    for idx = ascSel
        set(posMarker,"XData",ascH(idx),"YData",ascM(idx));
        addpoints(traceLine,ascH(idx),ascM(idx));
        [v,gifFrameIndex,frameSize] = emitFrame_(v,f,gifPath,gifFrameIndex,frameSize,options.AlsoWriteGif);
    end

    plot(axLeft,ascH,ascM,"Color",RECORDED_COLOR,"LineWidth",1.1,"HandleVisibility","off")
    clearpoints(traceLine)

    set(caption,"String",sprintf("curve %d of %d",k,nCurves))
    [v,gifFrameIndex,frameSize] = emitFrame_(v,f,gifPath,gifFrameIndex,frameSize,options.AlsoWriteGif);
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
    [v,gifFrameIndex,frameSize] = emitFrame_(v,f,gifPath,gifFrameIndex,frameSize,options.AlsoWriteGif);
end

hold(axLeft,"off")
hold(axRight,"off")

savePlot(f,baseName + "_finalFrame.pdf")

% Close explicitly (rather than waiting for the onCleanup at function
% exit) so the AVI is fully flushed to disk before ffmpeg tries to read
% it below.
close(v)

mp4Path = fullfile(plotsDir,baseName + ".mp4");
[ffmpegStatus,~] = system("ffmpeg -version");
if ffmpegStatus == 0
    transcodeCmd = sprintf('ffmpeg -y -loglevel error -i "%s" -c:v libx264 -pix_fmt yuv420p -movflags +faststart "%s"', ...
        aviPath,mp4Path);
    [transcodeStatus,transcodeMsg] = system(transcodeCmd);
    if transcodeStatus == 0 && isfile(mp4Path)
        delete(aviPath)
        videoPath = mp4Path;
    else
        warning("animateFORCIdealHysteron:TranscodeFailed", ...
            "ffmpeg transcode to MP4 failed; keeping the Motion JPEG AVI instead.\n%s",transcodeMsg);
        videoPath = aviPath;
    end
else
    videoPath = aviPath;
    fprintf("ffmpeg not found on the system path; delivering Motion JPEG AVI instead of MP4.\n");
end

fprintf("Wrote %s\n",videoPath);
if options.AlsoWriteGif
    fprintf("Wrote %s\n",gifPath);
end

end


function idx = subsampleIndices_(n,maxFrames)
%SUBSAMPLEINDICES_ Evenly spaced indices 1:n, capped at maxFrames, always including both ends.
idx = unique(round(linspace(1,n,min(maxFrames,n))));
end


function [v,frameIndex,frameSize] = emitFrame_(v,f,gifPath,frameIndex,frameSize,writeGif)
%EMITFRAME_ Capture the current figure into the video, and optionally the GIF.
%
%   Uses GETFRAME, which turned out not to be the problem: what actually
%   corrupted this animation's frames (stray diagonal lines, garbled
%   text) was recording with VideoWriter's default MPEG-4/H.264 profile
%   in this environment -- confirmed with an independent decoder
%   (ffmpeg), so it was a broken bitstream, not a display/capture
%   artifact. Once the caller records to Motion JPEG (see the top of
%   this file) instead, plain GETFRAME is fine and far cheaper than
%   round-tripping every frame through EXPORTGRAPHICS and a temp file.
%
%   GETFRAME's pixel size can still drift if the figure's content
%   bounding box changes (e.g. the colorbar/legend added partway
%   through this animation), and VideoWriter requires every frame to
%   match the first one it was given. FRAMESIZE is fixed from this
%   function's first call and every later frame is resized to match.
drawnow
frame = getframe(f);
img = frame.cdata;
if isempty(frameSize)
    frameSize = [size(img,1) size(img,2)];
elseif ~isequal([size(img,1) size(img,2)],frameSize)
    img = imresize(img,frameSize);
end
writeVideo(v,img);
if writeGif
    [indexed,map] = rgb2ind(img,256);
    if frameIndex == 0
        imwrite(indexed,map,gifPath,"gif","LoopCount",Inf,"DelayTime",1/max(v.FrameRate,1));
    else
        imwrite(indexed,map,gifPath,"gif","WriteMode","append","DelayTime",1/max(v.FrameRate,1));
    end
end
frameIndex = frameIndex + 1;
end


function closeIfOpen_(v)
%CLOSEIFOPEN_ Close a VideoWriter, tolerating one already closed explicitly.
%   This is the onCleanup safety net for an error mid-animation; the
%   normal path already closes v itself before attempting the ffmpeg
%   transcode, so a second close() here is expected to no-op.
try
    close(v)
catch
end
end
