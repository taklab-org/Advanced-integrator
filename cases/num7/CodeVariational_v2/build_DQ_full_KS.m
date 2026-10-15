function DQ_full = build_DQ_full_KS(mu_k_orig, Nop_coeff, a_bar, h, n_cheb, lambda_sup, L_sp)
% BUILD_DQ_FULL_KS  Build the full DQ operator for 1D KS (sine basis).
%
%   Uses BLOCK ordering: index = (s-1)*n_cheb + (n+1)
%
%   Builds DQ by computing N_op_lin_KS for each basis element column-by-column,
%   using the same padding convention as v2's F_lin_fin_dim_KS.

Ma = size(a_bar);
N_a = Ma - 1;

% Spatial mode sizes from mu_k_orig
if isvector(mu_k_orig)
    num_sp = numel(mu_k_orig);
else
    num_sp = prod(size(mu_k_orig));
end
total_dim = num_sp * n_cheb;

% K_sq computation
k_vals = (0:num_sp-1)';
j_max = size(Nop_coeff, 1) - 1;
K_q_vals = zeros(num_sp, 1);
for s = 1:num_sp
    K_q_vals(s) = (k_vals(s) * L_sp)^j_max;
end

% Pad a_bar to Ma+1 (matching v2's F_lin_fin_dim_KS convention)
a_bar_pad = zeros(Ma + 1);
for i = 1:length(Ma)
    otherdims{i} = 1:Ma(i);
end
a_bar_pad(otherdims{:}) = a_bar;

% Build DQ column by column
DQ_full = zeros(total_dim, total_dim);

% Part 1: lambda_sup * I
for idx = 1:total_dim
    DQ_full(idx, idx) = lambda_sup;
end

% Part 2: Nonlinear operator column-by-column
pad_sz = Ma + 1;
for s2 = 1:num_sp
    for ell2 = 0:n_cheb-1
        col_idx = (s2 - 1) * n_cheb + ell2 + 1;

        % Build unit vector e_{s2, ell2} in padded space
        e_pad = zeros(pad_sz);
        e_pad(s2, ell2 + 1) = 1;

        % Compute N_op_lin_KS(Nop_coeff, a_bar_pad, e_pad)
        nop_col = N_op_lin_KS(Nop_coeff, a_bar_pad, e_pad);

        % Apply K_sq and extract entries for each spatial mode
        for s1 = 1:num_sp
            K_q = K_q_vals(s1);
            if K_q == 0
                continue;
            end
            for ell1 = 0:n_cheb-1
                row_idx = (s1 - 1) * n_cheb + ell1 + 1;
                val = K_q * nop_col(s1, ell1 + 1);
                if abs(val) > 0
                    DQ_full(row_idx, col_idx) = DQ_full(row_idx, col_idx) + val;
                end
            end
        end
    end
end

end
