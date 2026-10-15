function [Psi, r, Y_0, Z_0, Z_1] = rig_1D_KS_variational_cheb(ba, m, Nop_coeff, h, params, nu, L)
% RIG_1D_KS_VARIATIONAL_CHEB  1D Kuramoto-Sivashinsky variational verification
%
%   Uses v2-compatible Newton solver (newton_F_KS) for sine-basis equation,
%   with v3 F-formulation bounds using build_DQ_full_KS_int (interval sine DQ).

n_cheb = size(ba, 1);
num_sp = m;

a_bar = ba.';

% Eigenvalues
c_dummy = zeros(m, n_cheb);
mu_k = mu_k_1DKS(size(c_dummy)-1, params, L);
mu_k_sp = mu_k(:, 1);

% =========================================================================
% Newton solve (v2-compatible, sine basis)
% =========================================================================
disp('=== Solving for fundamental matrix columns (KS Newton) ===');
Psi = zeros(num_sp * n_cheb, num_sp);

for i = 1:num_sp
    e_j = zeros(m, 1);
    e_j(i) = 1;
    [c_sol, ~] = newton_F_KS(mu_k, Nop_coeff, a_bar, zeros(m, n_cheb), h, e_j, L);
    Psi(:, i) = reshape(c_sol.', [], 1);
end

% =========================================================================
% Bounds (F-formulation, sine-basis interval DQ)
% =========================================================================
disp('=== Computing bounds context ===');
[~, ~, Z_1, Z_0, bounds_ctx] = compute_bounds_cheb(mu_k_sp, Nop_coeff, a_bar, h, nu, L, m, @build_DQ_full_KS_int, false);

if Z_1 >= 1
    error('Z_1 = %g >= 1. Cannot verify. Try smaller h or larger n_cheb.', Z_1);
end

Y_0 = zeros(1, num_sp);
for i = 1:num_sp
    beta = zeros(num_sp * n_cheb, 1);
    beta((i-1)*n_cheb + 1) = 1;
    [Y_i, ~, ~] = compute_Y_cheb(Psi(:, i), beta, bounds_ctx, nu, h);
    Y_0(i) = Y_i;
end

den = 1 - Z_1;
r = max(Y_0) / den;

disp('=== Results ===');
disp(['  Z_0     = ', num2str(Z_0)]);
disp(['  Z_inf   = ', num2str(bounds_ctx.Z_inf)]);
disp(['  Z_1     = ', num2str(Z_1)]);
disp(['  max(Y_0)= ', num2str(max(Y_0))]);
disp(['  r       = ', num2str(r)]);

end
