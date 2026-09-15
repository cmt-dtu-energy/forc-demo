function d = forcDistribution(Hr,H,M,smoothingFactor)
%FORCDISTRIBUTION Pike's mixed second derivative of a FORC family.
%
%   D = FORCDISTRIBUTION(HR,H,M,SMOOTHINGFACTOR) fits
%
%       M ~ a1 + a2*u + a3*u^2 + a4*v + a5*v^2 + a6*u*v
%
%   by least squares over the (2*SF+1)-by-(2*SF+1) block of grid points
%   around each point of the diagram, with u and v the offsets in Hr
%   and H, and reads the mixed second derivative rho = -1/2*a6 off the
%   fit. Noisy families cannot be differentiated twice by finite
%   differences alone, which is why this is a fit rather than a
%   difference of neighbouring points.
%
%   SF is both the smoothing width and the resolution limit: a larger
%   SF gives a cleaner but blunter distribution. rho is NaN wherever
%   the window does not fit entirely inside the measured wedge
%   H >= Hr (the diagonal and outer edges).
%
%   D is a plain struct with fields Hr, H, rho, Hc, Hu and
%   SmoothingFactor. Hc = (H-Hr)/2 and Hu = (H+Hr)/2 are returned as
%   matrices, one entry per (Hr,H) grid point, since that rotated grid
%   is curvilinear.
%
%   This is the one copy of the fit; FORCSimulation.distribution and
%   the plotFORC*.m scripts all call this instead of keeping their own
%   local copy.
%
%   See also FORCSimulation.

arguments (Input)
    Hr (1,:) double
    H (1,:) double
    M double
    smoothingFactor (1,1) double {mustBePositive, mustBeInteger}
end

sf = smoothingFactor;
nHr = numel(Hr);
nH = numel(H);
rho = NaN(nHr,nH);
for i = (1+sf):(nHr-sf)
    rowsIdx = (i-sf):(i+sf);
    for j = (1+sf):(nH-sf)
        colsIdx = (j-sf):(j+sf);

        block = M(rowsIdx,colsIdx);
        if any(isnan(block(:)))
            continue
        end

        [V,U] = meshgrid(H(colsIdx) - H(j),Hr(rowsIdx) - Hr(i));
        u = U(:);
        v = V(:);

        su = max(abs(u));
        sv = max(abs(v));
        if su == 0 || sv == 0
            continue
        end
        u = u/su;
        v = v/sv;

        A = [ones(numel(u),1), u, u.^2, v, v.^2, u.*v];
        coefficients = A\block(:);

        rho(i,j) = -0.5*coefficients(6)/(su*sv);
    end
end

[HH,HHr] = meshgrid(H,Hr);
d = struct("Hr",Hr,"H",H,"rho",rho, ...
    "Hc",(HH - HHr)/2,"Hu",(HH + HHr)/2, ...
    "SmoothingFactor",sf);

end
