function F = F_map(a,a0,n,N,L,params,h)
% a includes zero mode
% Input: size(a)  = [N1*N2*n, 1] : Vector form (One-sided Chebyshev*cosFourier*cosFourier, 1)
% Output: size(F) = [N1*N2*n, 1] : Vector form ((One-sided Chebyshev + initial defect)*cosFourier*cosFourier, 1)
%
% params = [lambda,sigma,epsilon];
% lambda  = params(1);
% sigma   = params(2);
% epsilon = params(3);
% Wave number: N
% N1 = N(1); N2 = N(2);
N1N2 = prod(N);
%
% Convert Vector form -> Matrix form: (cosFourier*cosFourier, One-sided Chebyshev)
a = reshape(a,[N1N2,n]);
%
% time derivative
da = (chebdiff_oneside(a.',h)).'; % Matrix form: (cosFourier*cosFourier, One-sided Chebyshev)
%
% Vector fields of SH 2D
F = RHS_OK_2D(a,n,N,L,params); % Matrix form: (cosFourier*cosFourier, One-sided Chebyshev)
%
% F1: defect of PDEs
F1 = da - F(:,1:n-1);
% F1 = -F(:,1:n-1);
% F1(2:end,:) = 2*F1(2:end,:); % back to Two-sided Chebyshev
%
% F2: defect of initial value
F2 = F2_map(a,a0);
% F2 = zeros(size(F2_map2(a,a0)));
% F: Matrix form: ((One-sided Chebyshev + initial defect), One-sided Chebyshev)
F = [F1, F2];
%
F = F(:);
