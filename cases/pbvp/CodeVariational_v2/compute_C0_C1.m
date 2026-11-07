function [C0_val, C1_val, rho_val] = compute_C0_C1(N, mu)
% COMPUTE_C0_C1  Compute the decay constants C_0(mu) and C_1(mu) from (24).
%
%   [C0_val, C1_val, rho_val] = compute_C0_C1(N, mu) computes:
%     C_0(mu): ||W_0 L_k^{-1} Lambda W_0^{-1}||_1 <= C_0(mu)/(2+mu)
%     C_1(mu): sup_n C_{1,n}(mu) where
%              1/(2-delta_{0,n}) sum_j |(W_0 L_k^{-1} Lambda)_{j,n}| <= C_{1,n}(mu)/(n+1+sqrt(n^2+mu^2))
%
%   Uses the finite truncation L_k^{(N)} with rigorous error control via rho_k.
%   Requires interval arithmetic (INTLAB).

if ~isa(mu, 'intval')
    mu = intval(mu);
end

% Step 1: Compute rho_k and L_k^{(N)} inverse
[rho_val, Linv_N] = compute_rho_k(N, mu);

if sup(rho_val) >= 1
    error('rho_k = %g >= 1: increase N for this mu value', sup(rho_val));
end

% Step 2: Build Lambda (conversion operator)
Lam = intval(conversion_op(N));

% Step 3: Compute W_0 * L_k^{(N)}^{-1} * Lambda * W_0^{-1}
% W_0 = diag(1, 2, 2, ..., 2), W_0^{-1} = diag(1, 1/2, 1/2, ..., 1/2)
w = 2 * intval(ones(N+1, 1));
w(1) = intval(1);
w_inv = intval(ones(N+1, 1)) / intval(2);
w_inv(1) = intval(1);

W0 = diag(w);
W0inv = diag(w_inv);

M = W0 * Linv_N * Lam * W0inv;

% Finite part: ||W_0 L_k^{(N)}^{-1} Lambda W_0^{-1}||_1
M_1norm = norm_1_mat(M);

% C_0^{(N)}: the norm from the finite truncation
% ||W_0 L_k^{-1} Lambda W_0^{-1}||_1 <= M_1norm / (1 - rho_val)
%   (accounting for truncation error via Neumann series)
C0_finite = M_1norm / (1 - rho_val);
C0_val = sup(C0_finite * (2 + mu));

% Step 4: Compute C_1(mu)
% C_{1,n}(mu) = (n+1+sqrt(n^2+mu^2)) / (2-delta_{0,n}) * sum_j |(W_0 L_k^{-1} Lambda)_{j,n}|
% We compute sup_{n=0,...,N} C_{1,n} over the FINITE truncation.
%
% NOTE: This C_1 is used only for n <= N2 (finite-to-tail bound Z_{<=N}).
%       Since N_for_C >= N2 (set in compute_bounds_cheb.m), the finite
%       columns 0,...,N2 are included in this supremum, so the bound is
%       safe for the current usage.
%       If C_1 were used as a UNIFORM bound over ALL n (including n > N),
%       additional verification would be required (see paper Lemma 2.5).

W0_Linv_Lam = W0 * Linv_N * Lam;

C1_max = intval(0);
for n = 0:N
    col_sum = sum(abs(W0_Linv_Lam(:, n+1)));
    if n == 0
        factor = 1;  % 1/(2 - delta_{0,0}) = 1/(2-1) = 1
    else
        factor = intval(1) / intval(2);  % 1/(2 - 0) = 1/2
    end
    alpha_n = intval(n) + 1 + sqrt(intval(n)^2 + mu^2);
    C1_n = alpha_n * factor * col_sum;
    C1_max = max(C1_max, C1_n);
end

% Account for truncation error
C1_val = sup(C1_max / (1 - rho_val));

end

function nrm = norm_1_mat(M)
% Compute the 1-norm (max column sum of absolute values) of a matrix
    n = size(M, 2);
    nrm = intval(0);
    for j = 1:n
        col_sum = sum(abs(M(:, j)));
        nrm = max(nrm, col_sum);
    end
end
