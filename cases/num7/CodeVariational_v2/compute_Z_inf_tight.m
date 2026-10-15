function Z_inf = compute_Z_inf_tight(mu_paper, C0_vals, h_iv, n_cheb, DQ_col_bound)
% COMPUTE_Z_INF_TIGHT  Tight Z_inf bound using tail-specific L_k^{-1} Lambda norm.
%
%   Z_inf = compute_Z_inf_tight(mu_paper, C0_vals, h_iv, n_cheb, DQ_col_bound)
%
%   Instead of bounding Z_inf = ||pi^{>N2} K pi^{>N2}|| using the FULL
%   operator norm ||L_k^{-1} Lambda|| = C_0/(2+mu) ~ O(1), this function
%   computes the TAIL-SPECIFIC norm:
%
%     ||pi^{>N2} L_k^{-1} Lambda||_{1,omega} ~ O(1/N2)
%
%   Key insight: the FULL norm C_0/(2+mu) is dominated by column 0
%   (where Chebyshev mode 0 maps strongly to modes 0,1,...). For the TAIL
%   (rows n > N2), the L_k^{-1} entries decay as 1/(2n), giving a column
%   sum ~ 1/N2 instead of ~ 1.

N2 = n_cheb - 1;
num_sp = length(mu_paper);

% Deduplicate: group spatial modes with the same mu value
mu_unique = unique_intval(mu_paper);
n_unique = length(mu_unique);

max_tail_norm = intval(0);

for u = 1:n_unique
    mu_s = mu_unique(u);

    % Find original index for C_0 lookup
    s_orig = find_mu_index(mu_paper, mu_s);
    C0_s = intval(C0_vals(s_orig));

    % Extended window: adapt to mu_s so that rho < 1
    N_ext = max(N2 + 20, ceil(2 * sup(abs(mu_s))) + 10);

    % =================================================================
    % Build L_k^{(N_ext)} and compute inverse (interval)
    % =================================================================
    L_ext = complete_op(N_ext, mu_s);
    Linv_ext = inv(L_ext);  % interval (N_ext+1) x (N_ext+1)

    % =================================================================
    % Compute tail norm using Lambda structure (avoid full matrix product)
    %
    % Lambda column j (1-based):
    %   row j+1: -1  (if j+1 <= N_ext+1)
    %   row j-1: +1  (if j >= 3)
    %
    % So M(:,j) = Linv*Lambda(:,j) = Linv(:,j-1) - Linv(:,j+1)  [j>=3]
    %           = -Linv(:,2)                                      [j=1]
    %           = -Linv(:,3)                                      [j=2]
    % =================================================================

    % Chebyshev weights: w_0 = 1, w_n = 2 for n >= 1
    w = 2 * intval(ones(N_ext+1, 1));
    w(1) = intval(1);

    % Tail row indices (1-based): N2+2 to N_ext+1
    tail_start = N2 + 2;  % 1-based index for Chebyshev mode N2+1 (0-based)

    tail_norm_s = intval(0);

    for j = 1:N_ext+1
        % Compute tail sum for column j of M = Linv * Lambda
        col_tail_sum = intval(0);

        if j == 1
            % M(:,1) = -Linv(:,2)
            for n = tail_start:N_ext+1
                col_tail_sum = col_tail_sum + w(n) * abs(Linv_ext(n, 2));
            end
        elseif j == 2
            % M(:,2) = -Linv(:,3)
            % (row j-1=1 has +1 only if j>=3, so column 2 only has row 3: -1)
            for n = tail_start:N_ext+1
                col_tail_sum = col_tail_sum + w(n) * abs(Linv_ext(n, 3));
            end
        elseif j <= N_ext
            % M(:,j) = Linv(:,j-1) - Linv(:,j+1)
            for n = tail_start:N_ext+1
                val = Linv_ext(n, j-1) - Linv_ext(n, j+1);
                col_tail_sum = col_tail_sum + w(n) * abs(val);
            end
        else
            % j = N_ext+1: M(:,end) = Linv(:,end-1)  (Lambda last row: -v_{N-1})
            for n = tail_start:N_ext+1
                col_tail_sum = col_tail_sum + w(n) * abs(Linv_ext(n, N_ext));
            end
        end

        tail_norm_s = max(tail_norm_s, col_tail_sum / w(j));
    end

    % =================================================================
    % Correction for truncation at N_ext
    % =================================================================
    % rho correction for finite truncation error
    [rho_ext, ~] = compute_rho_k(N_ext, mu_s);
    if sup(rho_ext) >= 1
        error('rho_ext >= 1 at N_ext=%d, mu=%.4g. Increase N_ext.', N_ext, sup(mu_s));
    end
    tail_norm_corrected = tail_norm_s / (1 - rho_ext);

    % Columns beyond N_ext: bounded by C_0(mu) / alpha_{N_ext+1}
    alpha_Nex = intval(N_ext+1) + 1 + sqrt(intval(N_ext+1)^2 + mu_s^2);
    tail_norm_beyond = C0_s / alpha_Nex;

    tail_norm_total = max(tail_norm_corrected, tail_norm_beyond);
    max_tail_norm = max(max_tail_norm, tail_norm_total);
end

Z_inf = h_iv/2 * max_tail_norm * DQ_col_bound;
Z_inf = sup(Z_inf);
Z_inf = intval(Z_inf);

end

% =========================================================================
% Helper: find unique intval values (by midpoint)
% =========================================================================
function mu_u = unique_intval(mu_paper)
    mids = mid(mu_paper(:));
    [~, ia] = unique(round(mids, 10));
    mu_u = mu_paper(ia);
end

% =========================================================================
% Helper: find index of matching mu value
% =========================================================================
function idx = find_mu_index(mu_paper, mu_target)
    mids = mid(mu_paper(:));
    target = mid(mu_target);
    [~, idx] = min(abs(mids - target));
end
