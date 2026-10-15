function [r, Y, Z_1, Z_0, details] = compute_bounds_cheb(mu_k_sp, Nop_coeff, a_bar, h, nu, L_sp, m, build_DQ_int_func, use_shift)
% COMPUTE_BOUNDS_CHEB  Compute radii-polynomial bounds (F-formulation).
%
%   [r, Y, Z_1, Z_0, details] = compute_bounds_cheb(mu_k_sp, Nop_coeff, a_bar, h, nu, L_sp, m, build_DQ_int_func, use_shift)
%
%   Implements bounds from Cadiot & Lessard, using the F-formulation:
%     F(U) = U + (h/2)*L^{-1}*Lambda*DQ*U - L^{-1}*beta   [eq (13)]
%     DF(U) = I + K,  K = (h/2)*L^{-1}*Lambda*DQ            [eq (21)]
%     A0 ~ inv(I + K)                                        [approx inverse of DF]
%
%   Simplified for:
%   - Finite-dimensional Fourier (Y_1 = Y_2 = 0)
%   - Linear variational problem (Z_2 = 0)
%   - Z_inf computed (NOT zero: Chebyshev tail-to-tail term remains)
%
%   Norm: ||U||_{1,omega} = sum_k omega_k |W_0 U_k|_1   [eq (16)]
%     W_0 = diag(1, 2, 2, ...) (no nu^n on Chebyshev)
%     omega_k = prod(2^{k_d>0}) per spatial dimension
%
%   Inputs:
%     mu_k_sp           : eigenvalues per spatial mode (flat vector)
%     Nop_coeff         : nonlinear operator coefficients
%     a_bar             : approximate PDE solution [m1, m2, ..., n_cheb]
%     h                 : time step
%     nu                : (unused in bounds, kept for interface compatibility)
%     L_sp              : spatial frequency scaling
%     m                 : spatial mode sizes
%     build_DQ_int_func : (optional) interval DQ builder, default @build_DQ_full_int
%     use_shift         : (optional, default true) If false, skip lambda_sup shift
%                         (uses mu_k = -(h/2)*lambda_k, DQ = N'(a_bar) only)

if ~exist('build_DQ_int_func', 'var') || isempty(build_DQ_int_func)
    build_DQ_int_func = @build_DQ_full_int;
end
if ~exist('use_shift', 'var') || isempty(use_shift)
    use_shift = true;
end

% Convert to intervals
h_iv = intval(h);
L_sp_iv = intval(L_sp);

% Problem dimensions
Ma = size(a_bar);
n_cheb = Ma(end);
N2 = n_cheb - 1;
num_sp = length(mu_k_sp);
total_dim = num_sp * n_cheb;

% Flatten eigenvalues
lambda_k = mu_k_sp(:);
lambda_k_iv = intval(lambda_k);

if use_shift
    lambda_sup = sup(max(lambda_k_iv));
    lambda_sup_iv = intval(lambda_sup);
    % Paper's mu_k: mu_k = (h/2)*(lambda_sup - lambda_k) >= 0
    mu_paper = h_iv/2 * (lambda_sup_iv - lambda_k_iv);
else
    lambda_sup_iv = intval(0);
    % Direct formulation: mu_k = -(h/2)*lambda_k (can be negative)
    mu_paper = -h_iv/2 * lambda_k_iv;
end

% Spatial mode sizes
Mc_sp = m;
if isscalar(Mc_sp)
    mu_k_sp_md = mu_k_sp(:);
else
    mu_k_sp_md = reshape(mu_k_sp, Mc_sp);
end

% =========================================================================
% Step 1: Build operators (interval arithmetic)
% =========================================================================
disp('Building operators (interval)...');

[L_full_iv, Lam_full_iv, ~] = build_Lk_full(intval(mu_k_sp), h_iv, n_cheb, use_shift);

% DQ in rigorous interval arithmetic
DQ_full_iv = build_DQ_int_func(intval(mu_k_sp_md), intval(Nop_coeff), intval(a_bar), h_iv, n_cheb, lambda_sup_iv, L_sp_iv);

% =========================================================================
% Step 2: Build K = (h/2)*L^{-1}*Lambda*DQ   [eq (21)]
% =========================================================================
disp('Computing K = (h/2)*L^{-1}*Lambda*DQ...');
Linv_iv = inv(L_full_iv);
clear L_full_iv;  % memory: only Linv needed downstream
K_iv = h_iv/2 * Linv_iv * Lam_full_iv * DQ_full_iv;

% DF = I + K
DF_iv = K_iv;
for ii = 1:total_dim
    DF_iv(ii,ii) = DF_iv(ii,ii) + intval(1);
end
clear K_iv;  % memory: K not used downstream

% =========================================================================
% Step 3: Compute A_0 = inv(I + K)  [approximate inverse of DF]
% =========================================================================
disp('Computing approximate inverse A_0...');
A0 = inv(mid(DF_iv));
A0_iv = intval(A0);

% =========================================================================
% Step 4: Compute Z_0 = ||I - A_0 * (I+K)||_{1,omega}   [Lemma 2.3]
% =========================================================================
disp('Computing Z_0...');
B = -A0_iv * DF_iv;
for ii = 1:total_dim
    B(ii,ii) = B(ii,ii) + intval(1);
end
clear DF_iv;  % memory: not used downstream
Mc = [Mc_sp, n_cheb];
Z_0 = compute_weighted_opnorm(B, Mc);
Z_0 = sup(Z_0);
clear B;  % memory: only Z_0 scalar needed
disp(['  Z_0 = ', num2str(Z_0)]);

% =========================================================================
% Step 5: Compute C_0, C_1 for each spatial mode
% =========================================================================
disp('Computing C_0, C_1 bounds...');
max_mu = sup(max(abs(mu_paper)));
N_for_C = max(N2, max(40, ceil(2 * max_mu)));
C0_vals = zeros(num_sp, 1);
C1_vals = zeros(num_sp, 1);
rho_vals = zeros(num_sp, 1);

for s = 1:num_sp
    mu_s = mu_paper(s);
    [C0_s, C1_s, rho_s] = compute_C0_C1(N_for_C, mu_s);
    C0_vals(s) = C0_s;
    C1_vals(s) = C1_s;
    rho_vals(s) = sup(rho_s);
end

if any(rho_vals >= 1)
    error('Some rho_k >= 1. Increase N_for_C.');
end
disp(['  max rho_k = ', num2str(max(rho_vals))]);

% =========================================================================
% Step 6: Compute DQ_col_bound from interval DQ matrix (rigorous)
% =========================================================================
disp('Computing ||DQ||_{1,omega} from interval matrix...');
DQ_col_bound = compute_weighted_opnorm(DQ_full_iv, Mc);
disp(['  ||DQ||_{1,omega} = ', num2str(sup(DQ_col_bound))]);

% =========================================================================
% Step 7: Compute Z_inf (Chebyshev tail-to-tail)
% =========================================================================
disp('Computing Z_inf (tight tail-specific bound)...');
Z_inf = compute_Z_inf_tight(mu_paper, C0_vals, h_iv, n_cheb, DQ_col_bound);
disp(['  Z_inf = ', num2str(sup(Z_inf))]);

% =========================================================================
% Step 8: Compute D_{A_0} diagonal (for Z_{A_0} bound)
% =========================================================================
disp('Computing D_{A_0}...');
A0_Linv_Lam = A0_iv * Linv_iv * Lam_full_iv;
clear A0_iv;  % memory: only A0 (double) saved in details

omega = build_omega_vec(Mc);

DA0_diag = intval(zeros(total_dim, 1));
for col = 1:total_dim
    % Same positive weighted sum, evaluated by rigorous interval array ops.
    col_sum = sum((omega / omega(col)) .* abs(A0_Linv_Lam(:, col)));
    s = ceil(col / n_cheb);
    n = mod(col - 1, n_cheb);
    mu_s = mu_paper(s);
    alpha_n = intval(n) + 1 + sqrt(intval(n)^2 + mu_s^2);
    DA0_diag(col) = col_sum * alpha_n;
end
clear A0_Linv_Lam;  % memory: not used downstream

% =========================================================================
% Step 9: Compute Z_{<=N} and Z_{A_0}
% =========================================================================
disp('Computing Z_{<=N} and Z_{A_0}...');

Z_leN = compute_Z_leN(mu_paper, C1_vals, h_iv, n_cheb, DQ_col_bound);
disp(['  Z_{<=N} = ', num2str(sup(Z_leN))]);

Z_A0 = compute_Z_A0(mu_paper, h_iv, n_cheb, DA0_diag, DQ_col_bound);
disp(['  Z_{A_0} = ', num2str(sup(Z_A0))]);

% =========================================================================
% Step 10: Compute Z_1   [eq (28)]
% =========================================================================
Z_0_iv = intval(Z_0);
% Paper eq (28): uses Z_inf^2, NOT Z_0^2
Z_1_term1 = Z_inf^2 + (1 + Z_inf)*Z_leN*Z_A0 + Z_A0*Z_inf^2;
Z_1_term2 = Z_0_iv*(1 + Z_leN*(1+Z_inf)) + Z_leN*Z_A0*(1 + Z_inf);
Z_1 = sup(max(Z_1_term1, Z_1_term2));
disp(['  Z_1 = ', num2str(Z_1)]);

if Z_1 >= 1
    warning('Z_1 = %g >= 1. Proof fails.', Z_1);
    r = inf;
    Y = inf;
    details = struct('Z_0', Z_0, 'Z_1', Z_1, 'Z_inf', sup(Z_inf), ...
                     'Z_leN', sup(Z_leN), 'Z_A0', sup(Z_A0));
    return;
end

% =========================================================================
% Step 11: Store context for Y computation
% =========================================================================
details = struct();
details.Z_0 = Z_0;
details.Z_1 = Z_1;
details.Z_inf = sup(Z_inf);
details.Z_leN = sup(Z_leN);
details.Z_A0 = sup(Z_A0);
details.C0_vals = C0_vals;
details.C1_vals = C1_vals;
details.rho_vals = rho_vals;
details.mu_paper = mu_paper;
details.A0 = A0;
details.Linv = Linv_iv;
details.Lam_full = Lam_full_iv;
details.DQ_full = DQ_full_iv;
details.Mc = Mc;
% Extra context for Y_3 direct evaluation
details.Nop_coeff = intval(Nop_coeff);
details.a_bar = intval(a_bar);
details.L_sp = L_sp_iv;
details.lambda_sup = lambda_sup_iv;
details.build_DQ_int_func = build_DQ_int_func;
details.n_fourier = Ma(1:end-1);  % spatial size of a_bar (may differ from m)

Y = nan;
r = nan;
disp('Bounds context computed successfully.');

end

% =========================================================================
% Helper functions
% =========================================================================

function nrm = compute_weighted_opnorm(B, Mc)
% Weighted operator norm ||B||_{1,omega} with W_0 weights (no nu^n)
    n_cheb = Mc(end);
    num_sp = prod(Mc(1:end-1));
    total_dim = num_sp * n_cheb;
    omega = build_omega_vec(Mc);

    nrm = intval(0);
    for col = 1:total_dim
        % INTLAB sum encloses the entire column; no unrounded double sum.
        weighted_col_sum = sum(omega .* abs(B(:, col)));
        nrm = max(nrm, weighted_col_sum / omega(col));
    end
end

function omega = build_omega_vec(Mc)
% Build omega weight vector: W_0 = {1,2,2,...} for Chebyshev,
% {1,2,2,...} per spatial dimension for cosine Fourier modes.
% No nu^n weighting on Chebyshev index.  [eq (16)]
    n_cheb = Mc(end);
    Mc_sp = Mc(1:end-1);
    num_sp = prod(Mc_sp);
    total_dim = num_sp * n_cheb;
    dim_sp = length(Mc_sp);

    sp_sz = Mc_sp;
    if isscalar(sp_sz), sp_sz = [sp_sz, 1]; end
    sp_subs = cell(1, max(dim_sp, 1));
    [sp_subs{:}] = ind2sub(sp_sz, (1:num_sp)');

    % Spatial weights: factor of 2 per dimension for k_d > 0
    w_sp = intval(ones(num_sp, 1));
    for d = 1:dim_sp
        for s = 1:num_sp
            if sp_subs{d}(s) > 1
                w_sp(s) = w_sp(s) * 2;
            end
        end
    end

    % Chebyshev weights: W_0 = {1, 2, 2, ...}  (NO nu^n)
    w_cheb = 2 * intval(ones(n_cheb, 1));
    w_cheb(1) = intval(1);

    omega = intval(zeros(total_dim, 1));
    for s = 1:num_sp
        for n = 1:n_cheb
            idx = (s-1)*n_cheb + n;
            omega(idx) = w_sp(s) * w_cheb(n);
        end
    end
end

function Z_inf = compute_Z_inf(mu_paper, C0_vals, h_iv, DQ_col_bound)
% Z_inf: Chebyshev tail-to-tail bound ||pi^{>N2} K pi^{>N2}||
% Bounded by C_0(mu_k)/(2+mu_k) * (h/2)*DQ_col_bound
    num_sp = length(mu_paper);
    Z_inf = intval(0);
    for s = 1:num_sp
        mu_s = mu_paper(s);
        C0_s = intval(C0_vals(s));
        tail_norm = C0_s / (2 + mu_s);
        contrib = h_iv/2 * tail_norm * DQ_col_bound;
        Z_inf = max(Z_inf, contrib);
    end
end

function Z_leN = compute_Z_leN(mu_paper, C1_vals, h_iv, n_cheb, DQ_col_bound)
% Z_{<=N}: tail-to-finite bound ||pi^{>N2} K pi^{<=N}||
% Bounded by C_1(mu_k)/alpha_{N2} * (h/2)*DQ_col_bound
    num_sp = length(mu_paper);
    N2 = n_cheb - 1;
    Z_leN = intval(0);
    for s = 1:num_sp
        mu_s = mu_paper(s);
        C1_s = intval(C1_vals(s));
        alpha_N2 = intval(N2) + 1 + sqrt(intval(N2)^2 + mu_s^2);
        tail_factor = C1_s / alpha_N2;
        contrib = h_iv/2 * tail_factor * DQ_col_bound;
        Z_leN = max(Z_leN, contrib);
    end
end

function Z_A0 = compute_Z_A0(mu_paper, h_iv, n_cheb, DA0_diag, DQ_col_bound)
% Z_{A_0}: A_0-weighted tail bound
    num_sp = length(mu_paper);
    N2 = n_cheb - 1;
    Z_A0 = intval(0);
    for s = 1:num_sp
        mu_s = mu_paper(s);
        for n = 0:N2
            idx = (s-1)*n_cheb + n + 1;
            alpha_n = intval(n) + 1 + sqrt(intval(n)^2 + mu_s^2);
            tail_decay = 1 / (intval(N2+1) + 1 + sqrt(intval(N2+1)^2 + mu_s^2));
            contrib = h_iv/2 * DA0_diag(idx) / alpha_n * tail_decay * DQ_col_bound;
            Z_A0 = max(Z_A0, contrib);
        end
    end
end

function nrm = norm_1_W0(a)
% ||a||_{1,W0}: sum over all elements with W_0 = {1,2,2,...} Chebyshev weights
    Ma = size(a);
    n_cheb = Ma(end);
    w = 2 * intval(ones(n_cheb, 1));
    w(1) = intval(1);
    nrm = intval(0);
    a_flat = reshape(a, [], n_cheb);
    for n = 1:n_cheb
        nrm = nrm + w(n) * sum(abs(a_flat(:, n)));
    end
end
