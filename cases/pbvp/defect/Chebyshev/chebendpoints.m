function epts = chebendpoints(a)% return function values of Chebyshev at endpoints
% a: Two-sided Chebyshev coefficients
%
% NOTE :: available for interval inputs, e.g., chebvals(intval(a),0,I))
%
w = (-1).^(0:length(a)-1);
epts = [w*a; sum(a)]; % fvals at -1, 1 or a, b
