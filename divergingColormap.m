function map = divergingColormap(n)
%DIVERGINGCOLORMAP A blue-white-red map for quantities with a meaningful zero.
%
%   MAP = divergingColormap(N) returns an N-by-3 colormap running from
%   blue through white to red.
%
%   A FORC distribution changes sign, and the sign is the whole point: a
%   negative lobe is a different statement about the sample than a weak
%   positive one. A sequential map such as parula hides that, since it
%   gives zero no distinguished colour. Used with a symmetric clim, white
%   here sits exactly at rho = 0.
%
%   See also createFigure, savePlot.

arguments (Input)
    n (1,1) double {mustBePositive, mustBeInteger} = 256
end

arguments (Output)
    map (:,3) double
end

low  = [0.02 0.19 0.38];    % deep blue
mid  = [1.00 1.00 1.00];    % white, at the centre of the scale
high = [0.40 0.00 0.12];    % deep red

half = floor(n/2);
lower = interp1([0;1], [low; mid], linspace(0,1,half)');
upper = interp1([0;1], [mid; high], linspace(0,1,n-half)');
map = [lower; upper];

end
