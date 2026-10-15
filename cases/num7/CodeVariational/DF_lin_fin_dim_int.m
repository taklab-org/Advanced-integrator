function [DF] = DF_lin_fin_dim_int(mu_k,Nop_coeff,a_bar,h,c,L,lambda)
% Block construction for interval Jacobian (from V2's DF_lin_fin_dim_int_v2).
% Replaces the O(numel_c^2) element-wise double loop with block construction.

if ~exist('lambda','var')
    lambda = [];
end

N = size(a_bar)-1;
P = size(Nop_coeff,2)-2;
Mc = size(c);
Ma = size(a_bar);
dim = length(Mc);
dim_sp = dim - 1;
n_cheb = Mc(end);
num_sp = prod(Mc(1:end-1));
numel_c = numel(mid(c));

% --- Preamble ---
mu_k_temp = mu_k;
mu_k_full = intval(zeros(Ma));
for i = 1:length(N)
    otherdims_test{i} = 1:N(i)+1;
    otherdims_c{i} = 1:Mc(i);
end
mu_k_full(otherdims_c{:}) = mu_k_temp;

a_bar_ext = intval(zeros((P+1)*N+1));
a_bar_ext(otherdims_test{:}) = a_bar;
N_ext = size(a_bar_ext)-1;
unit = intval(zeros(size(a_bar_ext)));
unit(1) = 1;

nop = N_op_lin_int(Nop_coeff, a_bar_ext, unit);
M_nop = size(nop)-1;

% K_sq computation (intval)
s_arr = arrayfun(@(kk) 0:kk, M_nop(1:end-1), 'UniformOutput', false);
k_cells = cell(1, dim_sp);
[k_cells{:}] = ndgrid(s_arr{:});
K_sq_arr = intval(zeros(size(k_cells{1})));
for ell = 1:dim_sp
    K_sq_arr = K_sq_arr + (k_cells{ell}*L(ell)).^(size(Nop_coeff,1)-1);
end
N_rep = ones(size(M_nop));
N_rep(end) = M_nop(end)+1;
K_sq_arr = repmat(K_sq_arr, N_rep);

nop_flat = reshape(nop, [numel(inf(a_bar_ext)), 1]);

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

% Per spatial mode: K_q (intval), mu_k spatial base index
K_q_vals = intval(zeros(num_sp, 1));
sp_mu_base = zeros(num_sp, 1);
for s = 1:num_sp
    K_q_vals(s) = K_sq_arr(1 + sum(k_sp(s,:) .* ext_strides));
    sp_mu_base(s) = 1 + sum(k_sp(s,:) .* N_strides);
end

% Lambda condition
skip_zero_sp = ~isempty(lambda);
is_zero_sp = all(k_sp == 0, 2);

% Chebyshev index matrices (double)
[E1, E2] = ndgrid(0:n_cheb-1, 0:n_cheb-1);

% Precompute template for Branches 1-2 (exact integers, single intval conversion)
block_12 = zeros(n_cheb, n_cheb);
block_12(1,1) = 1;
if n_cheb > 1
    block_12(1, 2:end) = 2*(-1).^(1:n_cheb-1);
end
for ell1 = 1:n_cheb-1
    block_12(ell1+1, ell1+1) = 2*ell1;
end
block_12_iv = intval(block_12);

% [Opt1] Pre-allocate intval zero templates (avoid intval constructor in loops)
block_5_zero = intval(zeros(n_cheb, n_cheb));
nop_zero = intval(zeros(n_cheb));
iv_zero = intval(0);

% [Opt1] Precompute sym_dif for all spatial modes (depends only on s2, not s1)
sp_signs_all = cell(num_sp, 1);
for s2 = 1:num_sp
    sp_signs_all{s2} = sym_dif(k_sp(s2,:));
end

% Build DF block by block
DF = intval(zeros(numel_c, numel_c));

for s1 = 1:num_sp
    idx_r = s1 + (0:n_cheb-1)*num_sp;
    active = ~skip_zero_sp || ~is_zero_sp(s1);

    % Branches 1-4: diagonal block
    block_diag = block_12_iv;
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
        for s2 = 1:num_sp
            idx_c = s2 + (0:n_cheb-1)*num_sp;
            sp_signs = sp_signs_all{s2};
            block_5 = block_5_zero;
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
                    nop_p = nop_zero;
                    nop_m = nop_zero;
                    lin_p = sp_base + cheb_stride*idx_p;
                    lin_m = sp_base + cheb_stride*idx_m;
                    nop_p(valid) = nop_flat(lin_p(valid));
                    nop_m(valid) = nop_flat(lin_m(valid));
                    contrib = h/2 * K_q * (nop_p - nop_m);
                    contrib(1,:) = iv_zero;
                    if sigma3 == -1
                        contrib(:,1) = iv_zero;
                    end
                    block_5 = block_5 + contrib;
                end
            end
            DF(idx_r, idx_c) = DF(idx_r, idx_c) + block_5;
        end
    end
end

DF = reshape(DF, [Mc, Mc]);
end
