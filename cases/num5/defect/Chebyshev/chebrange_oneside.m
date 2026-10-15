function [range] = chebrange_oneside(a,I)% return mag (max of abs value) of Chebyshev
% a: One-sided Chebyshev coefficients
% I: domain of the Chebyshev
%
% NOTE :: available for interval inputs, e.g., chebmag(intval(a),I))
%
arguments
    a
    I = [-1,1];
end
a = [a(1,:); 2*a(2:end,:)];
[maxima,rmax,minima,rmin] = chebminmax(a,I);
range = hull(minima,maxima);
% %
% if exist('intval','file') && isintval(a(1))
%     mag = intval(zeros(ndim,1));
%     r = intval(zeros(ndim,1));
%     for i=1:ndim
%         [m1,r1] = chebmag_main(a(:,i),I);
%         mag(i) = m1; r(i) = r1;
%     end
% else
%     mag = zeros(ndim,1);
%     r = zeros(ndim,1);
%     for i=1:ndim
%         [m1,r1] = chebmag_main(a(:,i),I);
%         mag(i) = m1; r(i) = r1;
%     end
% end
% end
% %
% function [maxima,rmax] = chebmag_main(a,I)
% % fvals at endpoints
% epts = chebendpoints(a);
% %
% % extrema of f
% r = chebextrema(a,I); rcand = [r;I'];
% if exist('intval','file') && isintval(a(1))
%     % find maximal value of Chebyshev using interval arithmetic
%     candidate = [chebvals(a,r,I);epts];
%     [maxima, indmax] = max(mag(candidate));
%     rmax = rcand(indmax);
% else
%     candidate = [chebvals(a,r,I);epts];
%     [maxima, indmax] = max(abs(candidate));
%     rmax = rcand(indmax);
% end
% end
