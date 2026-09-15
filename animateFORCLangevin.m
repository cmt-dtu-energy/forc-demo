function animateFORCLangevin(outputFile,options)
%ANIMATEFORCLANGEVIN Animate an equilibrium Langevin FORC family.

arguments (Input)
    outputFile (1,1) string = ""
    options.OutputFile (1,1) string = "langevinAnimation.mp4"
    options.FrameRate (1,1) double {mustBePositive} = 15
    options.MaxCurves (1,1) double {mustBePositive, mustBeInteger} = 12
    options.SmoothingFactor (1,1) double {mustBePositive, mustBeInteger} = 2
    options.DescentFrames (1,1) double {mustBePositive, mustBeInteger} = 8
    options.AscentFrames (1,1) double {mustBePositive, mustBeInteger} = 12
    options.HoldFramesAtEnd (1,1) double {mustBeNonnegative, mustBeInteger} = 40
    options.AlsoWriteGif (1,1) logical = false
end

if outputFile == ""
    outputFile = fullfile(fileparts(mfilename("fullpath")), ...
        "forc_demo_langevin_output.json");
end

raw = jsondecode(fileread(outputFile));
Hr = raw.diagram.Hr(:)';
H = raw.diagram.H(:)';
M = raw.diagram.M;
Hmajor = raw.majorBranch.H(:)';
Mmajor = raw.majorBranch.M(:)';
d = forcDistribution(Hr,H,M,options.SmoothingFactor);
peak = max(abs(d.rho),[],'all','omitnan');
if isempty(peak) || ~isfinite(peak) || peak == 0
    peak = 1;
end

nCurves = min(options.MaxCurves,numel(Hr));
rows = unique(round(linspace(1,numel(Hr),nCurves)));
plotsDir = fullfile(fileparts(mfilename("fullpath")),"plots");
if ~isfolder(plotsDir)
    mkdir(plotsDir)
end
videoPath = fullfile(plotsDir,options.OutputFile);
gifPath = fullfile(plotsDir,options.OutputFile + ".gif");

v = VideoWriter(videoPath,"MPEG-4");
v.FrameRate = options.FrameRate;
open(v);
closeVideo = onCleanup(@() close(v));

f = figure;
f.Position(3:4) = [1500 700];
tl = tiledlayout(f,1,2,"TileSpacing","compact","Padding","compact");
axLeft = nexttile(tl);
hold(axLeft,"on"); grid(axLeft,"on"); axis(axLeft,"square")
xlim(axLeft,[min(Hmajor) max(Hmajor)]);
ylim(axLeft,[min(Mmajor,[],'all') max(Mmajor,[],'all')]);
xlabel(axLeft,"Applied field, H"); ylabel(axLeft,"Magnetization, M");
title(axLeft,"Langevin field sweep and reversal curves");
plot(axLeft,Hmajor,Mmajor,"k","LineWidth",2,"DisplayName","Major descending branch");
marker = plot(axLeft,NaN,NaN,"o","MarkerSize",9,"MarkerFaceColor",[0.85 0.1 0.1], ...
    "MarkerEdgeColor","k","HandleVisibility","off");
trace = animatedline(axLeft,"Color",[0.85 0.1 0.1],"LineWidth",2,"HandleVisibility","off");
legend(axLeft,"Location","northwest");

axRight = nexttile(tl);
hold(axRight,"on"); grid(axRight,"on"); axis(axRight,"square")
xlim(axRight,[min(d.Hc,[],'all') max(d.Hc,[],'all')]);
ylim(axRight,[min(d.Hu,[],'all') max(d.Hu,[],'all')]);
xlabel(axRight,"Coercivity coordinate, H_c");
ylabel(axRight,"Interaction coordinate, H_u");
title(axRight,"FORC distribution (revealed after the family)");
caption = text(axRight,mean(xlim(axRight)),mean(ylim(axRight)), ...
    sprintf("curve 0 of %d",numel(rows)),"HorizontalAlignment","center");

frameIndex = 0;
for k = 1:numel(rows)
    row = rows(k);
    if k == 1
        startH = Hmajor(1);
    else
        startH = Hr(rows(k-1));
    end
    [~,idxStart] = min(abs(Hmajor-startH));
    [~,idxEnd] = min(abs(Hmajor-Hr(row)));
    descH = Hmajor(idxStart:idxEnd);
    descM = Mmajor(idxStart:idxEnd);
    for idx = subsample_(numel(descH),options.DescentFrames)
        set(marker,"XData",descH(idx),"YData",descM(idx));
        addpoints(trace,descH(idx),descM(idx));
        [v,frameIndex] = frame_(v,f,gifPath,frameIndex,options.AlsoWriteGif);
    end
    cols = find(~isnan(M(row,:)));
    ascH = H(cols);
    ascM = M(row,cols);
    for idx = subsample_(numel(ascH),options.AscentFrames)
        set(marker,"XData",ascH(idx),"YData",ascM(idx));
        addpoints(trace,ascH(idx),ascM(idx));
        [v,frameIndex] = frame_(v,f,gifPath,frameIndex,options.AlsoWriteGif);
    end
    plot(axLeft,ascH,ascM,"Color",[0.35 0.35 0.35],"HandleVisibility","off");
    clearpoints(trace);
    set(caption,"String",sprintf("curve %d of %d",k,numel(rows)));
    [v,frameIndex] = frame_(v,f,gifPath,frameIndex,options.AlsoWriteGif);
end

delete(caption);
contourf(axRight,d.Hc,d.Hu,d.rho,24,"LineColor","none","HandleVisibility","off");
colormap(axRight,divergingColormap);
clim(axRight,[-peak peak]);
yline(axRight,0,"LineStyle",":","Color","k","HandleVisibility","off");
colorbar(axRight);
title(axRight,sprintf("Equilibrium Langevin FORC distribution, SF = %d", ...
    options.SmoothingFactor));
for i = 1:options.HoldFramesAtEnd
    [v,frameIndex] = frame_(v,f,gifPath,frameIndex,options.AlsoWriteGif);
end

[~,base] = fileparts(options.OutputFile);
savePlot(f,base + "_finalFrame.pdf");
close(f);
fprintf("Wrote %s\n",videoPath);

end

function indices = subsample_(n,maxFrames)
indices = unique(round(linspace(1,n,min(maxFrames,n))));
end

function [video,frameIndex] = frame_(video,figureHandle,gifPath,frameIndex,writeGif)
drawnow;
frame = getframe(figureHandle);
writeVideo(video,frame);
if writeGif
    [indexed,map] = rgb2ind(frame.cdata,256);
    if frameIndex == 0
        imwrite(indexed,map,gifPath,"gif","LoopCount",Inf,"DelayTime",1/video.FrameRate);
    else
        imwrite(indexed,map,gifPath,"gif","WriteMode","append","DelayTime",1/video.FrameRate);
    end
end
frameIndex = frameIndex + 1;
end
