function [maxima,rmax,minima,rmin] = chebminmax(a,I)% return maxima/minima of Chebyshev & these coordinates
% a: Two-sided Chebyshev coefficients
% I: domain of the Chebyshev
%
% NOTE :: available for interval inputs, e.g., chebminmax(intval(a),I))
%
arguments
    a
    I = [-1,1];
end
% dimension of Chebyshev interpolations
ndim = size(a,2);
%
if exist('intval','file') && isintval(a(1))
    maxima = intval(zeros(ndim,1));
    rmax = intval(zeros(ndim,1));
    minima = intval(zeros(ndim,1));
    rmin = intval(zeros(ndim,1));
    for i=1:ndim
        [m1,r1,m2,r2] = chebminmax_main(a(:,i),I);
        maxima(i) = m1; rmax(i) = r1;
        minima(i) = m2; rmin(i) = r2;
    end
else
    maxima = zeros(ndim,1);
    rmax = zeros(ndim,1);
    minima = zeros(ndim,1);
    rmin = zeros(ndim,1);
    for i=1:ndim
        [m1,r1,m2,r2] = chebminmax_main(a(:,i),I);
        maxima(i) = m1; rmax(i) = r1;
        minima(i) = m2; rmin(i) = r2;
    end
end

end

function [maxima,rmax,minima,rmin] = chebminmax_main(a,I)
% fvals at endpoints
epts = chebendpoints(a);
%
% extrema of f
r = chebextrema(a,I); rcand = [r;I'];
if exist('intval','file') && isintval(a(1))
    % find minimal value of Chebyshev using interval arithmetic
    candidate = [chebvals(a,r,I);epts];
    [~, indmin] = min(inf(candidate));
    minima = candidate(indmin);
    rmin = rcand(indmin);
    % find maximal value of Chebyshev using interval arithmetic
    [~, indmax] = max(sup(candidate));
    maxima = candidate(indmax);
    rmax = rcand(indmax);
else
    candidate = [chebvals(a,r,I);epts];
    [minima, indmin] = min(candidate);
    rmin = rcand(indmin);
    [maxima, indmax] = max(candidate);
    rmax = rcand(indmax);
end
end
