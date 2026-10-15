function ChebCoeffs = chebcoeffs(f,M,I)% return Two-sided Chebyshev
% f: given function handle
% M: size of Chebyshev (maximum wave # is M-1)
% I: domain of the function f
arguments
    f
    M = 101
    I = [-1,1]
end
a = I(1); b = I(2);
% n = M-1;
cpts  = chebpts(M, a, b);
fvals = f(cpts);
FourierCoeffs = (fft([fvals;flipud(fvals(2:end-1,:))]));
ChebCoeffs = FourierCoeffs(1:M,:)/(M-1);
ChebCoeffs(1,:) = ChebCoeffs(1,:)/2;
ChebCoeffs(end,:) = ChebCoeffs(end,:)/2;

% Post-process:
if ( isreal(fvals) )
    % Real-valued case:
    ChebCoeffs = real(ChebCoeffs);
elseif ( isreal(1i*fvals) )
    % Imaginary-valued case:
    ChebCoeffs = 1i*imag(ChebCoeffs);
end
end
