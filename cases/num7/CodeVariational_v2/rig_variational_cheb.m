function [Psi, r, Y_0, Z_0, Z_1] = rig_variational_cheb(ba, m, Nop_coeff, h, lambda, nu, L_sp, mu_k_func, build_DQ_func, build_DQ_int_func, use_shift)
% RIG_VARIATIONAL_CHEB  Rigorous verification for variational equations
%                        using the upgraded Chebyshev integration method.
%
%   [Psi, r, Y_0, Z_0, Z_1] = rig_variational_cheb(ba, m, Nop_coeff, h, ...
%       lambda, nu, L_sp, mu_k_func, build_DQ_func, build_DQ_int_func, use_shift)
%
%   Inputs:
%     ba              : approximate PDE solution [n_cheb, m1, m2, ...]
%     m               : spatial mode sizes [m1, m2, ...]
%     Nop_coeff       : nonlinear operator coefficients
%     h               : time step
%     lambda          : physical parameter (passed to mu_k_func)
%     nu              : weight parameter (kept for interface, unused in bounds)
%     L_sp            : spatial frequency scaling
%     mu_k_func       : function handle for eigenvalues
%     build_DQ_func   : (optional) double DQ builder, default @build_DQ_full
%     build_DQ_int_func: (optional) interval DQ builder, default @build_DQ_full_int
%     use_shift       : (optional, default true) If false, skip lambda_sup shift.
%                        For LINEAR variational equations (Z_2=0), no-shift
%                        is valid and reduces ||DQ|| by removing lambda_sup*I.

if ~exist('build_DQ_func', 'var') || isempty(build_DQ_func)
    build_DQ_func = @build_DQ_full;
end
if ~exist('build_DQ_int_func', 'var') || isempty(build_DQ_int_func)
    build_DQ_int_func = @build_DQ_full_int;
end
if ~exist('use_shift', 'var') || isempty(use_shift)
    use_shift = true;
end

% Reshape a_bar: from [n_cheb, m1, m2, ...] to [m1, m2, ..., n_cheb]
a_bar = permute(ba, [2:length(size(ba)), 1]);
n_cheb = size(ba, 1);
num_sp = prod(m);

% Compute eigenvalues
c_dummy = zeros([m, n_cheb]);
mu_k_full = mu_k_func(size(c_dummy)-1, lambda, L_sp);

idx_cell = cell(1, length(m));
for d = 1:length(m)
    idx_cell{d} = 1:m(d);
end
idx_cell{end+1} = 1;
mu_k_sp = mu_k_full(idx_cell{:});
mu_k_sp = reshape(mu_k_sp, [], 1);

if use_shift
    lambda_sup = max(mu_k_sp);
else
    lambda_sup = 0;
end

if isscalar(m)
    mu_k_sp_md = mu_k_sp;
else
    mu_k_sp_md = reshape(mu_k_sp, m);
end

% Build operators (double precision for Newton)
disp('=== Building operators (double) ===');
[L_full, Lam_full, ~] = build_Lk_full(mu_k_sp, h, n_cheb, use_shift);
DQ_full = build_DQ_func(mu_k_sp_md, Nop_coeff, a_bar, h, n_cheb, lambda_sup, L_sp);

% =========================================================================
% Precompute bounds context (F-formulation, interval arithmetic)
% =========================================================================
disp('=== Computing bounds context ===');
[~, ~, Z_1, Z_0, bounds_ctx] = compute_bounds_cheb(mu_k_sp, Nop_coeff, a_bar, h, nu, L_sp, m, build_DQ_int_func, use_shift);

if Z_1 >= 1
    error('Z_1 = %g >= 1. Cannot verify. Try smaller h or larger n_cheb.', Z_1);
end

% =========================================================================
% Solve for each column of the fundamental matrix
% =========================================================================
disp('=== Solving for fundamental matrix columns ===');
Psi = zeros(num_sp * n_cheb, num_sp);
Y_0 = zeros(1, num_sp);

for i = 1:num_sp
    beta = zeros(num_sp * n_cheb, 1);
    beta((i-1)*n_cheb + 1) = 1;

    c0 = zeros(num_sp * n_cheb, 1);
    [c_sol, ~] = newton_cheb(L_full, Lam_full, DQ_full, beta, h, c0);

    [Y_i, ~, ~] = compute_Y_cheb(c_sol, beta, bounds_ctx, nu, h);
    Y_0(i) = Y_i;

    Psi(:, i) = c_sol;

    if mod(i, max(1, floor(num_sp/10))) == 0
        disp(['  Column ', num2str(i), '/', num2str(num_sp), ...
               ', Y = ', num2str(Y_i)]);
    end
end

% =========================================================================
% Final bounds
% =========================================================================
den = 1 - Z_1;
r = max(Y_0) / den;

disp('=== Results ===');
disp(['  Z_0     = ', num2str(Z_0)]);
disp(['  Z_inf   = ', num2str(bounds_ctx.Z_inf)]);
disp(['  Z_1     = ', num2str(Z_1)]);
disp(['  max(Y_0)= ', num2str(max(Y_0))]);
disp(['  r       = ', num2str(r)]);
disp(['  1-Z_1   = ', num2str(den)]);

end
