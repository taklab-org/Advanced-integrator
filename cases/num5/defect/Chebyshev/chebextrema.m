function ext = chebextrema(a,I)% return coordinates of extrema for Chebyshev series
% a: Two-sided Chebyshev coefficients
% I: domain of the Chebyshev
%
% NOTE :: available for interval inputs
%
arguments
    a
    I = [-1,1];
end
ext = chebroots(chebdiff(a),I);
