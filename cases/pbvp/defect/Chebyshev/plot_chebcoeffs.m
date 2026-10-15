function plot_chebcoeffs(a)
% a: given (two-sided) chebyshev coefficients
[M,~] = size(a); % M: size of Chebyshev
% N = (m+1)/2;
k = 0:(M-1);
semilogy(k,abs(a),"LineWidth",1.6)
xlabel('$k$','interpreter','latex')
ylabel('$|a_{k,\ell}|$','interpreter', 'latex')
