function F = F_map_linear(a,ba,a0,n,N,L,params,h)
% Input: size(a)  = [N*n, 1] : Vector form (One-sided Chebyshev*cosFourier, 1)
% Output: size(F) = [N*n, 1] : Vector form ((One-sided Chebyshev + initial defect)*cosFourier, 1)
%
% Wave number: N
N1N2 = prod(N);
%
% Convert Vector form -> Matrix form: (cosFourier, One-sided Chebyshev)
a = reshape(a,[N1N2,n]);
%
% time derivative
da = (chebdiff_oneside(a.',h)).'; % Matrix form: (cosFourier, One-sided Chebyshev)
%
% Vector fields of SH 1D
F = RHS_OK_2D_linear(a,ba,n,N,L,params); % Matrix form: (cosFourier, One-sided Chebyshev)
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
