function L = complete_op(N, mu)
% BUILD_COMPLETE_OP  Build the finite truncation L_k^{(N)} of the operator L_k.
%
%   L = complete_op(N, mu) returns an (N+1)x(N+1) matrix representing the
%   Chebyshev integration operator L_k from the paper (eq. 11):
%
%     Row 0 (initial condition): [1, -2, 2, -2, 2, ...]
%     Row j (j>=1):              [..., mu, 2j, -mu, ...]
%
%   Here mu = (h/2)*(lambda_sup - lambda_k) >= 0.
%
%   Works with both double and intval inputs.

L = zeros(N+1, N+1);
if isa(mu, 'intval')
    L = intval(L);
end

% Row 0: initial condition  u_k(-1) = b_k
% T_n(-1) = (-1)^n, so u_k(-1) = U_{k,0} + 2*sum_{n>=1} (-1)^n U_{k,n}
L(1,1) = 1;
for j = 1:N
    L(1, j+1) = 2*(-1)^j;
end

% Rows 1..N: tridiagonal structure
% Row j: mu*c_{j-1} + 2j*c_j - mu*c_{j+1}
for j = 1:N
    L(j+1, j+1) = 2*j;       % diagonal: 2j
    L(j+1, j) = mu;           % sub-diagonal: mu
    if j < N
        L(j+1, j+2) = -mu;   % super-diagonal: -mu
    end
end

end
