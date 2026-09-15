function [f,ax] = createFigure(fs)
% CREATEFIGURE Create and configure a figure with a gridded axis ready to hold plots
arguments (Input)
    fs (1,1) double {mustBePositive,mustBeReal,mustBeFinite} % fontsize for all elements in the figure
end

arguments (Output)
    f (1,1) matlab.ui.Figure
    ax (1,1) matlab.graphics.axis.Axes
end

    f = figure;
    ax = gca;
    fontsize(f,fs,"points")
    hold(ax,"on")
    grid(ax,"on");
    axis square
end
