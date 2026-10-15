function x = chebpts(M, a, b)
% M: size of Chebyshev polynomials (M-1 is maximum order of Chebyshev)
arguments
    M; a=-1; b=1;
end
    tt = linspace(0, pi, M).';
    xi = cos(tt);
    x = (1.0 - xi).*a/2 + (1.0 + xi).*b/2;
end
