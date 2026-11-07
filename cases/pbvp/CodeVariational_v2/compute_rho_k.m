function [rho, Linv_N] = compute_rho_k(N, mu)
% COMPUTE_RHO_K  Compute rho_k from Lemma 3.1 to verify L_k invertibility.
%
%   [rho, Linv_N] = compute_rho_k(N, mu) computes:
%     rho_k = |mu_k|/2 * max{rho^(1), rho^(2), rho^(3)}
%
%   where rho^(1), rho^(2), rho^(3) are defined in Lemma 3.1.
%   If rho_k < 1, then L_k is boundedly invertible.
%
%   Also returns Linv_N = (L_k^{(N)})^{-1} as interval matrix.
%
%   Uses interval arithmetic for rigorous bounds.

L_N = complete_op(N, mu);
Linv_N = inv(L_N);

% First column of Linv_N (W_0 weighted)
col0 = Linv_N(:, 1);
w = 2 * intval(ones(N+1, 1));
w(1) = intval(1);
W0_col0_norm = sum(abs(col0) .* w);

% Column N+1 of Linv_N (last column, 0-indexed: column N)
colN = Linv_N(:, N+1);
W0_colN1_norm = sum(abs(colN) .* w);

% rho^(1): 1/(N+1) * (||W_0 * col(0)||_1 + 1)
rho1 = 1 / intval(N+1) * (W0_col0_norm + 1);

% rho^(2): ||W_0 * col(N+1)||_1 + 1/(N+2) * ||W_0 * col(0)||_1 + 1/(N+2)
rho2 = W0_colN1_norm + 1/intval(N+2) * W0_col0_norm + 1/intval(N+2);

% rho^(3): 2/((N+1)*(N+3)) * ||W_0 * col(0)||_1 + 1/(N+1) + 1/(N+3)
rho3 = 2 / (intval(N+1)*intval(N+3)) * W0_col0_norm + 1/intval(N+1) + 1/intval(N+3);

% rho_k = |mu|/2 * max{rho1, rho2, rho3}
rho = abs(mu) / 2 * max(max(rho1, rho2), rho3);

end
