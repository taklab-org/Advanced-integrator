function [DF] = DF_lin_fin_dim(mu_k,Nop_coeff,a_bar,h,c,L,lambda)
% Fused version: V2 block algorithm + New1 persistent cache.
% - First call: builds DF using O(num_sp^2) block construction (V2)
% - Subsequent calls with same parameters: returns cached result instantly (New1)
% DF does not depend on c values (only size(c)), so caching is valid.

if ~exist('lambda','var')
    lambda = [];
end

% --- Persistent cache (from New1) ---
persistent cache
if local_cache_hit(cache, mu_k, Nop_coeff, a_bar, h, size(c), L, lambda)
    DF = cache.DF;
    return
end

% --- V2 block construction algorithm ---
N = size(a_bar)-1;
P = size(Nop_coeff,2)-1;
Mc = size(c);
Ma = size(a_bar);
dim = length(Mc);
dim_sp = dim - 1;
n_cheb = Mc(end);
num_sp = prod(Mc(1:end-1));
numel_c = prod(Mc);

% Preamble
mu_k_temp = mu_k;
mu_k_full = zeros(Ma);
for i = 1:length(N)
    otherdims_test{i} = 1:N(i)+1;
    otherdims_c{i} = 1:Mc(i);
end
mu_k_full(otherdims_c{:}) = mu_k_temp;

a_bar_ext = zeros(P*N+1);
a_bar_ext(otherdims_test{:}) = a_bar;
N_ext = size(a_bar_ext)-1;
unit = zeros(size(a_bar_ext));
unit(1) = 1;

nop = N_op_lin(Nop_coeff, a_bar_ext, unit);
M_nop = size(nop)-1;

% K_sq computation
s_arr = arrayfun(@(kk) 0:kk, M_nop(1:end-1), 'UniformOutput', false);
k_cells = cell(1, dim_sp);
[k_cells{:}] = ndgrid(s_arr{:});
if (size(Nop_coeff,1)-1) > 0
    K_sq_arr = zeros(size(k_cells{1}));
    for ell = 1:dim_sp
        K_sq_arr = K_sq_arr + (k_cells{ell}*L(ell)).^(size(Nop_coeff,1)-1);
    end
else
    K_sq_arr = ones(size(k_cells{1}));
end
N_rep = ones(size(M_nop));
N_rep(end) = M_nop(end)+1;
K_sq_arr = repmat(K_sq_arr, N_rep);

nop_flat = reshape(nop, [], 1);

% Precompute spatial mode table
sp_subs = cell(1, dim_sp);
sp_sz = Mc(1:end-1); if isscalar(sp_sz), sp_sz = [sp_sz, 1]; end
[sp_subs{:}] = ind2sub(sp_sz, (1:num_sp)');
k_sp = zeros(num_sp, dim_sp);
for d = 1:dim_sp
    k_sp(:,d) = sp_subs{d} - 1;
end

% Strides for N_ext-sized array
ext_strides = ones(1, dim_sp);
for d = 2:dim_sp
    ext_strides(d) = prod(N_ext(1:d-1)+1);
end
cheb_stride = prod(N_ext(1:dim_sp)+1);

% Strides for N-sized array
N_strides = ones(1, dim_sp);
for d = 2:dim_sp
    N_strides(d) = prod(N(1:d-1)+1);
end
cheb_stride_N = prod(N(1:dim_sp)+1);

% Per spatial mode: K_q, mu_k base index
K_q_vals = zeros(num_sp, 1);
sp_mu_base = zeros(num_sp, 1);
for s = 1:num_sp
    K_q_vals(s) = K_sq_arr(1 + sum(k_sp(s,:) .* ext_strides));
    sp_mu_base(s) = 1 + sum(k_sp(s,:) .* N_strides);
end

% Lambda condition
skip_zero_sp = ~isempty(lambda);
is_zero_sp = all(k_sp == 0, 2);

% Chebyshev index matrices
[E1, E2] = ndgrid(0:n_cheb-1, 0:n_cheb-1);

% Build DF block by block
DF = zeros(numel_c, numel_c);

for s1 = 1:num_sp
    idx_r = s1 + (0:n_cheb-1)*num_sp;
    active = ~skip_zero_sp || ~is_zero_sp(s1);

    % Branches 1-4: diagonal block
    block_diag = zeros(n_cheb, n_cheb);
    block_diag(1,1) = 1;
    if n_cheb > 1
        block_diag(1, 2:end) = 2*(-1).^(1:n_cheb-1);
    end
    for ell1 = 1:n_cheb-1
        block_diag(ell1+1, ell1+1) = 2*ell1;
    end
    if active
        for ell1 = 1:n_cheb-1
            mu_val = mu_k_full(sp_mu_base(s1) + ell1*cheb_stride_N);
            block_diag(ell1+1, ell1) = block_diag(ell1+1, ell1) - h/2*mu_val;
            if ell1 <= n_cheb-2
                block_diag(ell1+1, ell1+2) = block_diag(ell1+1, ell1+2) + h/2*mu_val;
            end
        end
    end
    DF(idx_r, idx_r) = block_diag;

    % Branch 5: nonlinear
    if active
        K_q = K_q_vals(s1);
        if K_q ~= 0
            for s2 = 1:num_sp
                idx_c = s2 + (0:n_cheb-1)*num_sp;
                sp_signs = sym_dif(k_sp(s2,:));
                block_5 = zeros(n_cheb, n_cheb);
                for m_idx = 1:size(sp_signs, 1)
                    sp_shift = abs(k_sp(s1,:) - sp_signs(m_idx,:));
                    if any(sp_shift > N_ext(1:dim_sp))
                        continue;
                    end
                    sp_base = 1 + sum(sp_shift .* ext_strides);
                    for sigma3 = [+1, -1]
                        cheb_diff = E1 - sigma3*E2;
                        idx_p = abs(cheb_diff + 1);
                        idx_m = abs(cheb_diff - 1);
                        valid = (idx_p <= N_ext(end)) & (idx_m <= N_ext(end));
                        nop_p = zeros(n_cheb);
                        nop_m = zeros(n_cheb);
                        lin_p = sp_base + cheb_stride*idx_p;
                        lin_m = sp_base + cheb_stride*idx_m;
                        nop_p(valid) = nop_flat(lin_p(valid));
                        nop_m(valid) = nop_flat(lin_m(valid));
                        contrib = h/2 * K_q * (nop_p - nop_m);
                        contrib(1,:) = 0;
                        if sigma3 == -1
                            contrib(:,1) = 0;
                        end
                        block_5 = block_5 + contrib;
                    end
                end
                DF(idx_r, idx_c) = DF(idx_r, idx_c) + block_5;
            end
        end
    end
end

DF = reshape(DF, [Mc, Mc]);

% --- Store in cache ---
cache = struct();
cache.mu_k = mu_k;
cache.Nop_coeff = Nop_coeff;
cache.a_bar = a_bar;
cache.h = h;
cache.c_size = size(c);
cache.L = L;
cache.lambda = lambda;
cache.DF = DF;
end

function hit = local_cache_hit(cache, mu_k, Nop_coeff, a_bar, h, c_size, L, lambda)
if isempty(cache)
    hit = false;
    return
end
hit = isequal(cache.c_size, c_size) ...
    && isequal(cache.h, h) ...
    && isequal(cache.L, L) ...
    && isequal(cache.lambda, lambda) ...
    && isequal(cache.mu_k, mu_k) ...
    && isequal(cache.Nop_coeff, Nop_coeff) ...
    && isequal(cache.a_bar, a_bar);
end
