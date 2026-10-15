function [Y, Y_0, Y_3] = compute_Y_cheb(c, e_j, details, nu, h)
% COMPUTE_Y_CHEB  Compute the Y bound (F-formulation).
%
%   [Y, Y_0, Y_3] = compute_Y_cheb(c, e_j, details, nu, h)
%
%   F-formulation:
%     F(c) = c + (h/2)*L^{-1}*Lambda*DQ*c - L^{-1}*beta   [eq (13)]
%     Y_0 = ||A_0 * pi^{<=N2} F(c)||_{1,omega}
%     Y_3 = ||pi^{>N2} L^{-1}*(rhs)||_{1,omega}   (Chebyshev tail)
%
%   rhs = (h/2)*Lambda*DQ*c - beta has FINITE support (up to P*N2+1)
%   because U_bar and c have finite Chebyshev support.
%   Therefore Y_3 is EVALUATED (not estimated) by solving on an extended
%   Chebyshev window and directly summing the tail entries.
%
%   For finite Fourier, linear problem (Y_1=Y_2=0):
%     Y = Y_0*(1 + Z_{<=N}*(1+Z_inf)) + Y_3*(1+Z_inf)*(Z_{A_0}+1)   [eq (32)]

h_iv = intval(h);
c_iv = intval(c);
e_j_iv = intval(e_j);

Mc = details.Mc;
n_cheb = Mc(end);
N2 = n_cheb - 1;
Mc_sp = Mc(1:end-1);
num_sp = prod(Mc_sp);
total_dim = num_sp * n_cheb;

mu_paper = details.mu_paper;
Linv_iv = details.Linv;

% =========================================================================
% Y_0: ||A_0 * pi^{<=N2} F(c)||_{1,omega}
% F(c) = c + (h/2)*L^{-1}*Lambda*DQ*c - L^{-1}*beta   [eq (13)]
% =========================================================================

F_iv = c_iv + h_iv/2 * Linv_iv * details.Lam_full * details.DQ_full * c_iv ...
       - Linv_iv * e_j_iv;
A0_F = intval(details.A0) * F_iv;

omega = build_omega_vec_local(Mc);
Y_0 = intval(0);
for i = 1:total_dim
    Y_0 = Y_0 + abs(A0_F(i)) * omega(i);
end
Y_0 = sup(Y_0);

% =========================================================================
% Y_3: DIRECT EVALUATION of pi^{>N2} L^{-1} * rhs
%
% rhs = (h/2)*Lambda*DQ*c - beta has finite support (up to P*N2+1)
% where P = max polynomial degree in Nop_coeff.
%
% Method:
%   1. Build rhs on extended Chebyshev window [0,...,n_ext-1]
%   2. For each spatial mode k, solve L_k^{(N_solve)} u = rhs_k
%   3. EVALUATE sum_{n=N2+1}^{N_solve} w_n |u_n| (interval arithmetic)
%   4. Bound tiny remainder beyond N_solve using truncation identity
% =========================================================================

Y_3 = compute_Y3_direct(c_iv, e_j_iv, details, h_iv);
Y_3 = sup(Y_3);

% =========================================================================
% Combine: Y = Y_0*(1 + Z_{<=N}*(1+Z_inf)) + Y_3*(1+Z_inf)*(Z_{A_0}+1)
% [eq (32) with Y_1=Y_2=0]
% =========================================================================
Z_leN = intval(details.Z_leN);
Z_A0 = intval(details.Z_A0);
Z_inf = intval(details.Z_inf);

Y = sup(intval(Y_0) * (1 + Z_leN*(1+Z_inf)) + intval(Y_3) * (1+Z_inf) * (Z_A0 + 1));

end

% =========================================================================
% Y_3 direct evaluation
% =========================================================================
function Y_3 = compute_Y3_direct(c_iv, beta_iv, details, h_iv)

    Mc = details.Mc;
    n_cheb = Mc(end);
    N2 = n_cheb - 1;
    Mc_sp = Mc(1:end-1);
    num_sp = prod(Mc_sp);
    mu_paper = details.mu_paper;

    % Nonlinearity degree P (max column index - 1 in Nop_coeff)
    P = size(details.Nop_coeff, 2) - 1;

    % Extended window: rhs support is at most P*N2 + 1 (from convolution + Lambda shift)
    n_ext = max(n_cheb, P * N2 + 2);

    % --- Build extended rhs = (h/2)*Lambda*DQ*c - beta ---
    rhs_ext = build_extended_rhs(c_iv, beta_iv, details, h_iv, n_ext);

    % Convert to per-mode blocks: [n_ext, num_sp]
    rhs_blocks = tensor_to_blocks(rhs_ext, Mc_sp, n_ext);

    % Spatial weights
    w_sp = build_spatial_weights(Mc_sp);

    % --- For each spatial mode: solve and evaluate tail ---
    Y_3 = intval(0);
    for s = 1:num_sp
        mu_s = mu_paper(s);
        rhs_s = rhs_blocks(:, s);

        tail_s = eval_tail_direct(rhs_s, mu_s, N2);
        Y_3 = Y_3 + w_sp(s) * tail_s;
    end
end

% =========================================================================
function rhs_ext = build_extended_rhs(c_iv, beta_iv, details, h_iv, n_ext)
% Build rhs = (h/2)*Lambda*(lambda_sup*c + K^j*N'(a_bar)*c) - beta
% on extended Chebyshev window of size n_ext.
%
% When n_fourier (spatial size of a_bar) differs from Mc_sp (= m),
% we zero-pad c to n_fourier for the convolution, then truncate back to m.

    Mc = details.Mc;
    n_cheb = Mc(end);
    Mc_sp = Mc(1:end-1);
    dim_sp = length(Mc_sp);

    % Spatial size of a_bar (may be larger than Mc_sp)
    if isfield(details, 'n_fourier')
        Nf_sp = details.n_fourier;
    else
        Nf_sp = Mc_sp;
    end

    % Reshape c and beta from block vectors to tensors [m1,...,md, n_cheb]
    c_tensor = blockvec_to_tensor(c_iv, Mc_sp, n_cheb);
    beta_tensor = blockvec_to_tensor(beta_iv, Mc_sp, n_cheb);

    % Zero-pad beta to extended Chebyshev window [m1,...,md, n_ext]
    ext_sz_m = [Mc_sp, n_ext];
    beta_ext = intval(zeros(ext_sz_m));
    low_idx_m = repmat({':'}, 1, length(ext_sz_m));
    low_idx_m{end} = 1:n_cheb;
    beta_ext(low_idx_m{:}) = beta_tensor;

    % For convolution: use n_fourier spatial size
    ext_sz_full = [Nf_sp, n_ext];
    c_ext_full = intval(zeros(ext_sz_full));
    a_ext_full = intval(zeros(ext_sz_full));

    % Embed c into full spatial window (zero-pad modes m+1..N)
    embed_idx = repmat({':'}, 1, length(ext_sz_full));
    for d = 1:dim_sp
        embed_idx{d} = 1:Mc_sp(d);
    end
    embed_idx{end} = 1:n_cheb;
    c_ext_full(embed_idx{:}) = c_tensor;

    % Fill a_bar into full spatial window
    low_idx_full = repmat({':'}, 1, length(ext_sz_full));
    low_idx_full{end} = 1:n_cheb;
    a_ext_full(low_idx_full{:}) = details.a_bar;

    % Compute N'(a_bar)*c on extended window (full spatial size)
    if isa(details.build_DQ_int_func, 'function_handle') && ...
       strcmp(func2str(details.build_DQ_int_func), 'build_DQ_full_KS_int')
        nop_ext_full = N_op_lin_KS_int(details.Nop_coeff, a_ext_full, c_ext_full);
    else
        nop_ext_full = N_op_lin_int(details.Nop_coeff, a_ext_full, c_ext_full);
    end

    % Truncate convolution result to m spatial modes
    trunc_idx = repmat({':'}, 1, length(ext_sz_full));
    for d = 1:dim_sp
        trunc_idx{d} = 1:Mc_sp(d);
    end
    nop_ext = nop_ext_full(trunc_idx{:});

    % Spatial multiplier K^j (differential operator powers) — on m modes
    K_mult = build_K_multiplier(ext_sz_m, details.Nop_coeff, details.L_sp);

    % c_ext on m modes for lambda_sup term
    c_ext = intval(zeros(ext_sz_m));
    c_ext(low_idx_m{:}) = c_tensor;

    % DQ*c = lambda_sup*c + K^j * N'(a_bar)*c
    dqc_ext = details.lambda_sup * c_ext + K_mult .* nop_ext;

    % Apply Lambda (conversion operator) in Chebyshev direction
    lam_dqc = apply_Lambda_cheb(dqc_ext);

    % rhs = (h/2)*Lambda*DQ*c - beta
    rhs_ext = h_iv / 2 * lam_dqc - beta_ext;
end

% =========================================================================
function tail_val = eval_tail_direct(rhs_s, mu_s, N2)
% Evaluate ||pi^{>N2} L_k^{-1} * rhs_s||_{W0} directly.
%
% Uses the FULL rhs (not just rhs_{>N2}) to correctly handle the case
% when c is an approximate Newton solution with nonzero residual.
%
% Previous version zeroed rhs_{<=N2} assuming pi^{>N2}(L_k^{-1}*rhs_{<=N2})=0,
% which holds only at the EXACT Newton solution. For a floating-point
% approximate solution with residual ~5e-11, the dropped term
% pi^{>N2}(L_k^{-1}*residual_{<=N2}) can be comparable to Y_3 itself.

    n_rhs = length(rhs_s);

    % Find support of full rhs
    support = find(sup(abs(rhs_s)) > 0, 1, 'last');
    if isempty(support)
        tail_val = intval(0);
        return;
    end

    % Choose N_solve large enough
    Nsolve = max([support + 10, N2 + 20, 50, ceil(2 * sup(abs(mu_s))) + 10]);

    % Build L_k^{(Nsolve)} and solve L*u = rhs (full rhs)
    L_N = complete_op(Nsolve, mu_s);
    Linv_N = inv(L_N);

    rhs_pad = intval(zeros(Nsolve + 1, 1));
    rhs_pad(1:min(n_rhs, Nsolve+1)) = rhs_s(1:min(n_rhs, Nsolve+1));
    u_N = Linv_N * rhs_pad;

    % Direct evaluation: sum |u_n| * w_n for n = N2+1,...,Nsolve
    tail_start = N2 + 2;  % 1-based
    if tail_start <= Nsolve + 1
        tail_eval = sum(2 * abs(u_N(tail_start:end)));
    else
        tail_eval = intval(0);
    end

    % Truncation remainder beyond Nsolve
    u_last = abs(u_N(Nsolve + 1));
    rho_tail = abs(mu_s) / (2 * intval(Nsolve + 2));

    if sup(rho_tail) >= 1
        error('rho_tail >= 1 at Nsolve=%d, mu=%.4g. Increase n_cheb.', Nsolve, sup(mu_s));
    end

    u_next = abs(mu_s) * u_last / (2 * intval(Nsolve + 1));
    tail_remainder = 2 * u_next / (1 - rho_tail);

    tail_val = tail_eval + tail_remainder;
end

% =========================================================================
% Helper: block vector → spatial tensor
% =========================================================================
function tensor = blockvec_to_tensor(vec, Mc_sp, n_cheb)
    num_sp = prod(Mc_sp);
    % vec is [num_sp * n_cheb, 1], stored as [cheb_block_1; cheb_block_2; ...]
    % Reshape to [n_cheb, num_sp] then to [n_cheb, m1, m2, ...]
    tmp = reshape(vec, [n_cheb, Mc_sp]);
    % Permute to [m1, m2, ..., n_cheb]
    nd = length(Mc_sp) + 1;
    tensor = permute(tmp, [2:nd, 1]);
end

% =========================================================================
% Helper: spatial tensor → per-mode blocks
% =========================================================================
function blocks = tensor_to_blocks(tensor, Mc_sp, n_ext)
    % tensor is [m1, m2, ..., n_ext]
    % Output: [n_ext, num_sp] where each column is one spatial mode
    nd = length(Mc_sp) + 1;
    % Permute to [n_ext, m1, m2, ...]
    perm = permute(tensor, [nd, 1:nd-1]);
    num_sp = prod(Mc_sp);
    blocks = reshape(perm, n_ext, num_sp);
end

% =========================================================================
% Helper: apply Lambda (conversion operator) along last dim
% =========================================================================
function Lv = apply_Lambda_cheb(v)
% (Lambda*v)_0 = 0
% (Lambda*v)_n = -v_{n-1} + v_{n+1}  for n = 1,...,N-2
% (Lambda*v)_N = -v_{N-1}             (truncated)
    sz = size(v);
    n_ext = sz(end);
    Lv = intval(zeros(sz));

    idx = repmat({':'}, 1, length(sz));

    for n = 2:n_ext-1
        cur = idx; cur{end} = n;
        prev = idx; prev{end} = n - 1;
        nxt = idx; nxt{end} = n + 1;
        Lv(cur{:}) = -v(prev{:}) + v(nxt{:});
    end

    % Last entry
    if n_ext >= 2
        cur = idx; cur{end} = n_ext;
        prev = idx; prev{end} = n_ext - 1;
        Lv(cur{:}) = -v(prev{:});
    end
end

% =========================================================================
% Helper: build spatial multiplier K^j
% =========================================================================
function K_mult = build_K_multiplier(ext_sz, Nop_coeff, L_sp)
    dim_sp = length(ext_sz) - 1;
    j_max = size(Nop_coeff, 1) - 1;

    if j_max == 0
        sp_shape = ext_sz(1:end-1);
        if isscalar(sp_shape), sp_shape = [sp_shape, 1]; end
        K_sp = intval(ones(sp_shape));
    else
        s_arr = arrayfun(@(kk) 0:kk-1, ext_sz(1:end-1), 'UniformOutput', false);
        k_cells = cell(1, dim_sp);
        [k_cells{:}] = ndgrid(s_arr{:});
        K_sp = intval(zeros(size(k_cells{1})));
        for ell = 1:dim_sp
            K_sp = K_sp + (intval(k_cells{ell}) * L_sp(ell)).^j_max;
        end
    end

    rep = ones(1, length(ext_sz));
    rep(end) = ext_sz(end);
    K_mult = repmat(K_sp, rep);
end

% =========================================================================
% Helper: spatial weights omega_k
% =========================================================================
function w_sp = build_spatial_weights(Mc_sp)
    num_sp = prod(Mc_sp);
    sp_sz = Mc_sp;
    if isscalar(sp_sz), sp_sz = [sp_sz, 1]; end
    sp_subs = cell(1, max(length(Mc_sp), 1));
    [sp_subs{:}] = ind2sub(sp_sz, (1:num_sp)');
    w_sp = intval(ones(num_sp, 1));
    for d = 1:length(Mc_sp)
        for s = 1:num_sp
            if sp_subs{d}(s) > 1
                w_sp(s) = w_sp(s) * 2;
            end
        end
    end
end

% =========================================================================
% Helper: build omega vector (for Y_0)
% =========================================================================
function omega = build_omega_vec_local(Mc)
    n_cheb = Mc(end);
    Mc_sp = Mc(1:end-1);
    num_sp = prod(Mc_sp);
    total_dim = num_sp * n_cheb;
    dim_sp = length(Mc_sp);

    sp_sz = Mc_sp;
    if isscalar(sp_sz), sp_sz = [sp_sz, 1]; end
    sp_subs = cell(1, max(dim_sp, 1));
    [sp_subs{:}] = ind2sub(sp_sz, (1:num_sp)');

    w_sp = intval(ones(num_sp, 1));
    for d = 1:dim_sp
        for s = 1:num_sp
            if sp_subs{d}(s) > 1
                w_sp(s) = w_sp(s) * 2;
            end
        end
    end

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
