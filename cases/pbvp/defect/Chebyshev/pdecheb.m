function cc = pdecheb(f, I, Wavenum, tol, Nmax)% return Two-sided Chebyshev
% f: ode function (Fourier coefficients) by ode45, ode78, ode89 etc.
% I: domain of the function f
% tol: tolerance of Chebyshev coefficients
% Nmax: maximum # of Chebyshev coefficients
arguments
    f
    I = [-1,1];
    Wavenum = (size(f.y,1)+1)/2 - 1; % 2N-1
    tol = eps;
    Nmax = 10000;
end
%
devalFun = @(x) deval(f, x).';
%
i = 3; N = (size(f.y,1)+1)/2; % 2N-1
%
while true
    test_size = 2^i+1;
    schbc = chebcoeffs(devalFun,test_size,I);
    chopsize = standardChop(schbc(:,N-Wavenum:N+Wavenum),tol);
    if test_size > chopsize || (2^i+1 > Nmax)
        break
    end
%     if all(abs(schbc(end-2:end,:)) < tol,'all') || (2^i+1 > Nmax)
%         break
%     end
    i = i + 1;
end
% chopsize
cc = schbc(1:chopsize,:);
%
end
