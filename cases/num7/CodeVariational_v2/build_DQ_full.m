function DQ_full = build_DQ_full(mu_k_orig, Nop_coeff, a_bar, h, n_cheb, lambda_sup, L_sp)
% BUILD_DQ_FULL  Build the full DQ operator (linearization of Q around a_bar).
%
%   Uses BLOCK ordering: index = (s-1)*n_cheb + (n+1)
%   where s = spatial mode (1..num_sp), n = Chebyshev mode (0..N2).
%
%   DQ * c represents: lambda_sup * c + N_op_linearized(a_bar, c)
%   (the "Q" operator from the paper, without h/2 or Lambda prefactors)

Ma = size(a_bar);
N_a = Ma - 1;
P = size(Nop_coeff, 2) - 1;
dim = length(Ma);
dim_sp = dim - 1;

% Spatial mode sizes from mu_k_orig (not a_bar, which may be larger)
if isvector(mu_k_orig)
    Mc_sp = numel(mu_k_orig);
else
    Mc_sp = size(mu_k_orig);
end
num_sp = prod(Mc_sp);
total_dim = num_sp * n_cheb;

% Extend a_bar for convolution
otherdims_test = cell(1, length(Ma));
for i = 1:length(Ma)
    otherdims_test{i} = 1:Ma(i);
end

a_bar_ext_sz = P * N_a + 1;
a_bar_ext = zeros(a_bar_ext_sz);
a_bar_ext(otherdims_test{:}) = a_bar;
N_ext = size(a_bar_ext) - 1;
unit = zeros(size(a_bar_ext));
unit(1) = 1;
nop = N_op_lin(Nop_coeff, a_bar_ext, unit);
M_nop = size(nop) - 1;

% K_sq computation
s_arr = arrayfun(@(kk) 0:kk, M_nop(1:end-1), 'UniformOutput', false);
k_cells = cell(1, dim_sp);
[k_cells{:}] = ndgrid(s_arr{:});
if (size(Nop_coeff,1)-1) > 0
    K_sq_arr = zeros(size(k_cells{1}));
    for ell = 1:dim_sp
        K_sq_arr = K_sq_arr + (k_cells{ell}*L_sp(ell)).^(size(Nop_coeff,1)-1);
    end
else
    K_sq_arr = ones(size(k_cells{1}));
end
N_rep = ones(size(M_nop));
N_rep(end) = M_nop(end)+1;
K_sq_arr = repmat(K_sq_arr, N_rep);
nop_flat = reshape(nop, [], 1);

% Spatial mode table
sp_sz = Mc_sp;
if isscalar(sp_sz), sp_sz = [sp_sz, 1]; end
sp_subs = cell(1, dim_sp);
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

K_q_vals = zeros(num_sp, 1);
for s = 1:num_sp
    K_q_vals(s) = K_sq_arr(1 + sum(k_sp(s,:) .* ext_strides));
end

% Chebyshev index matrices
[E1, E2] = ndgrid(0:n_cheb-1, 0:n_cheb-1);

% =========================================================================
% Build DQ matrix with BLOCK ordering
% Block (s1, s2): n_cheb x n_cheb sub-matrix
% Row index: (s1-1)*n_cheb + (ell1+1)
% Col index: (s2-1)*n_cheb + (ell2+1)
% =========================================================================
DQ_full = zeros(total_dim, total_dim);

% Part 1: lambda_sup * I (diagonal)
for s = 1:num_sp
    for n = 0:n_cheb-1
        idx = (s-1)*n_cheb + n + 1;
        DQ_full(idx, idx) = lambda_sup;
    end
end

% Part 2: Nonlinear operator linearization
% This is the Chebyshev convolution contribution from N_op_lin(a_bar, c)
% adapted from DF_lin_fin_dim block-5 structure.
%
% The N_op contribution to the equation is:
%   nop(a_bar, c)_{k, ell} = sum over spatial/Chebyshev convolution indices
%
% In the DQ formulation (without Lambda or h/2):
% DQ acts on c and produces the convolution term that will be hit by (h/2)*Lambda
for s1 = 1:num_sp
    K_q = K_q_vals(s1);
    if K_q == 0
        continue;
    end

    for s2 = 1:num_sp
        sp_signs = sym_dif(k_sp(s2,:));
        block_nop = zeros(n_cheb, n_cheb);

        for m_idx = 1:size(sp_signs, 1)
            sp_shift = abs(k_sp(s1,:) - sp_signs(m_idx,:));
            if any(sp_shift > N_ext(1:dim_sp))
                continue;
            end
            sp_base = 1 + sum(sp_shift .* ext_strides);

            for sigma3 = [+1, -1]
                cheb_diff = E1 - sigma3*E2;
                idx_conv = abs(cheb_diff);
                valid = (idx_conv <= N_ext(end));
                nop_vals = zeros(n_cheb);
                lin_idx = sp_base + cheb_stride * idx_conv;
                nop_vals(valid) = nop_flat(lin_idx(valid));

                contrib = K_q * nop_vals;

                if sigma3 == -1
                    contrib(:,1) = 0;
                end
                block_nop = block_nop + contrib;
            end
        end

        % Place block into DQ matrix
        idx_r = (s1-1)*n_cheb + (1:n_cheb);
        idx_c = (s2-1)*n_cheb + (1:n_cheb);
        DQ_full(idx_r, idx_c) = DQ_full(idx_r, idx_c) + block_nop;
    end
end

end
