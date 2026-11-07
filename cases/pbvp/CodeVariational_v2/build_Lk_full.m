function [L_full, Lam_full, lambda_sup] = build_Lk_full(mu_k_orig, h, n_cheb, use_shift)
% BUILD_LK_FULL  Build the full block-diagonal L and Lambda operators.
%
%   [L_full, Lam_full, lambda_sup] = build_Lk_full(mu_k_orig, h, n_cheb, use_shift)
%
%   Inputs:
%     mu_k_orig : array of eigenvalues lambda_k for each spatial mode
%     h         : time step
%     n_cheb    : number of Chebyshev modes
%     use_shift : (optional, default true) If true, use lambda_sup shift:
%                 mu_k = (h/2)*(lambda_sup - lambda_k) >= 0.
%                 If false, use direct eigenvalues:
%                 mu_k = -(h/2)*lambda_k (can be negative).
%
%   When use_shift=false, lambda_sup is returned as 0 (no shift applied).

if nargin < 4 || isempty(use_shift)
    use_shift = true;
end

lambda_k = mu_k_orig(:);  % flatten spatial modes
num_sp = length(lambda_k);
N = n_cheb - 1;  % Chebyshev truncation (0-indexed max)

if use_shift
    % Standard formulation: shift so all mu_k >= 0
    if isa(lambda_k, 'intval')
        lambda_sup = max(sup(lambda_k));
        lambda_sup = intval(lambda_sup);
    else
        lambda_sup = max(lambda_k);
    end
    mu_paper = (h/2) * (lambda_sup - lambda_k);  % >= 0 for all k
else
    % Direct formulation: no shift, mu_k = -(h/2)*lambda_k
    if isa(lambda_k, 'intval')
        lambda_sup = intval(0);
    else
        lambda_sup = 0;
    end
    mu_paper = -(h/2) * lambda_k;  % can be negative
end

% Build block-diagonal L and Lambda
total_dim = num_sp * n_cheb;
L_full = zeros(total_dim, total_dim);
Lam_full = zeros(total_dim, total_dim);

if isa(mu_paper, 'intval')
    L_full = intval(L_full);
    Lam_full = intval(Lam_full);
end

Lam_block = conversion_op(N);
if isa(mu_paper, 'intval')
    Lam_block = intval(Lam_block);
end

for s = 1:num_sp
    idx = (s-1)*n_cheb + (1:n_cheb);
    L_full(idx, idx) = complete_op(N, mu_paper(s));
    Lam_full(idx, idx) = Lam_block;
end

end
