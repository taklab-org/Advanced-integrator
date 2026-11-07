function Lam = conversion_op(N)
% CONVERSION_OP  Build the Chebyshev conversion/integration operator Lambda.
%
%   Lam = conversion_op(N) returns an (N+1)x(N+1) matrix representing
%   the operator Lambda from the paper (eq. 11):
%
%     (Lambda v)_0 = 0
%     (Lambda v)_n = -v_{n-1} + v_{n+1}  for n >= 1
%
%   With truncation at N: (Lambda v)_N = -v_{N-1}.

Lam = zeros(N+1, N+1);

for j = 1:N
    Lam(j+1, j) = -1;       % sub-diagonal: -1
    if j < N
        Lam(j+1, j+2) = 1;  % super-diagonal: +1
    end
end

end
