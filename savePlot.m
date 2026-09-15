function savePlot(f,fname)
% SAVEPLOT Save a figure into a local "plots" subfolder.
arguments (Input)
    f (1,1) matlab.ui.Figure % figure to save to file
    fname (1,1) string % filename under the plots directory to save the figure
end
    plotsDir = fullfile(fileparts(mfilename("fullpath")),"plots");
    if ~isfolder(plotsDir)
        mkdir(plotsDir);
    end
    exportgraphics(f, ...
        fullfile(plotsDir, fname), ...
        "Resolution",600);
end
