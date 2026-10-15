function T = tridiag_int(N, mu)
% TRIDIAG_INT  Build the tridiagonal part T_mu of L_k (without the initial condition row).
%
%   T = tridiag_int(N, mu) returns an NxN matrix:
%     T_mu = diag(mu, -1) + diag(2, 4, 6, ..., 2N) + diag(-mu, +1)
%
%   This is the core of L_k from rows 1..N (0-indexed), i.e. L_k without
%   the initial condition row.

T = zeros(N, N);
if isa(mu, 'intval')
    T = intval(T);
end

for j = 1:N
    T(j, j) = 2*j;          % diagonal: 2j
    if j > 1
        T(j, j-1) = mu;     % sub-diagonal: mu
    end
    if j < N
        T(j, j+1) = -mu;    % super-diagonal: -mu
    end
end

end
