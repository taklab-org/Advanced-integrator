function F = RHS_OK_2D_linear(a,ba,n,N,L,params)
% Input: size(a)  = [N1N2, n] : Matrix form (cosFourier*cosFourier, One-sided Chebyshev)
% Output: size(F) = [N1N2, n] : Matrix form (cosFourier*cosFourier, One-sided Chebyshev)
%
lambda  = params(1);
sigma   = params(2);
epsilon = params(3);
%
% Wave number: N
N1 = N(1); N2 = N(2);
L1 = L(1); L2 = L(2);
%
a = reshape(a,[N1,N2,n]);   % size(a)  = [N1,N2,n];
ba = reshape(ba,[N1,N2,n]); % size(ba) = [N1,N2,n];
%
% Derivative operator
% k = (0:N-1)';
[kx,ky] = ndgrid(0:N1-1,0:N2-1);
%
% Linear operator
bkL2 = (kx*L1).^2 + (ky*L2).^2;
mu_k = bkL2 .* (-epsilon^2*bkL2 + 1) - sigma;
%
% Cubic Convolution
% a3 = powerconvcos(a,3);
ba2a = cubicconvcos(ba,ba,a);
%
% Vector fields : One-sided Chebyshev
F = mu_k.*a - 3*bkL2.*ba2a;
%
F = reshape(F,N1*N2,n);
F(1,1) = F(1,1) + lambda;
end
