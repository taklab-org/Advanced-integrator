function F = F2_map(a,a0)
% a includes zero mode
% Input: size(a)  = [N1*N2, n] : Matrix form (cosFourier*cosFourier, One-sided Chebyshev)
%        size(a0) = [N1,N2]    : Matrix form (cosFourier,cosFourier)
% Output: size(F) = [N*N2, 1] : Vector form (cosFourier)

a(:, 2:end) = 2*a(:, 2:end); % One-sided --> Two-sided Chebyshev

n = size(a,2); % size of Chebyshev
ell = (0:n-1).';
v = (-1).^ell;

F = a*v - a0;
