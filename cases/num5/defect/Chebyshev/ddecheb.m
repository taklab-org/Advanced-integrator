function cc = ddecheb(f, I, tol, Nmax)% return Two-sided Chebyshev
% f: dde function by dde23, etc.
% I: domain of the function f
% tol: tolerance of Chebyshev coefficients
% Nmax: maximum # of Chebyshev coefficients
arguments
    f
    I = [-1,1];
    tol = eps;
    Nmax = 10000;
end
%
devalFun = @(x) deval(f, x).';
cc = cheb(devalFun,I,tol,Nmax);
end
