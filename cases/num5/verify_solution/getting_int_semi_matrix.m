function [int_semi_mat, int_semi_mat_q, ba2_sup, e_K_inf, e_inf_K, ba2_norm, finite] = getting_int_semi_matrix(ba, h, m, params, L, nu, directional_only)
% getting_int_semi_matrix - Rigorous finite integrated semigroup matrix majorant for 2D-OK
%
% Returns two m1m2 x m1m2 interval matrices:
%   int_semi_mat   - integrated semigroup (for Y-bound / defect propagation)
%   int_semi_mat_q - int_semi_mat * diag(bkL2) (for Z-bound / nonlinear Lipschitz)
% Optional directional_only=true skips scalar legacy outputs 4 and 5;
% the finite-flow data in output 7 and all directional bounds are unchanged.
%
% The OK equation has N(a)_k = -|k|^2_L * (a^3)_k, so the nonlinear
% perturbation carries a bkL2 = |k|^2_L factor. This requires separate
% semigroup evaluations for the defect bound (Y) and Lipschitz bound (Z).

sigma   = intval(params(2));
if nargin<8,directional_only=nargout<4;end
epsilon = intval(params(3));
m1 = m(1); m2 = m(2); m1m2 = prod(m);

[~, ba2] = powerconvcos(ba, 2);
size_ba2 = size(ba2);
n_ba2 = size_ba2(1); N1_ba2 = size_ba2(2); N2_ba2 = size_ba2(3);

ba2_sup = reshape(chebmag_oneside(reshape(ba2, n_ba2, N1_ba2*N2_ba2)), N1_ba2, N2_ba2);
ba2_norm = 3 * wnorm(ba2_sup, intval(nu));
% These legacy scalar diagnostics are not used by the directional coupling
% or the intermediate CAP. Avoid their costly scans when callers ignore them.
e_K_inf=intval(0);e_inf_K=intval(0);
if ~directional_only
    ba2_sup_tail=ba2_sup;ba2_sup_tail(1,1)=0;
    c_conv=3*wnorm(ba2_sup_tail,intval(nu));
    if nargout>=4,e_K_inf=min(3*getting_cm(ba2_sup,m,nu),c_conv);end
    if nargout>=5,e_inf_K=min(3*getting_cinf(ba2_sup,m,nu),c_conv);end
end

%% Build DN matrix using proper 2D cosine convolution
Chebyshev = size(ba2, 1);

% Mode indices for finite block (ndgrid order = column-major linearisation)
[kx_grid, ky_grid] = ndgrid(0:m1-1, 0:m2-1);
kx_vec = kx_grid(:);   % m1m2 x 1
ky_vec = ky_grid(:);

DN_cheb_cos_cos = intval(zeros(Chebyshev, m1m2, m1m2));
for j_lin = 1:m1m2          % source mode (l1,l2)
  l1 = kx_vec(j_lin); l2 = ky_vec(j_lin);
  for i_lin = 1:m1m2        % target mode (k1,k2)
    k1 = kx_vec(i_lin); k2 = ky_vec(i_lin);
    if l1 == 0 && l2 == 0
      DN_cheb_cos_cos(:,i_lin,j_lin) = 3 * ba2(:, k1+1, k2+1);
    elseif l2 == 0           % l1 > 0
      DN_cheb_cos_cos(:,i_lin,j_lin) = 3 * ( ...
        ba2(:, abs(k1-l1)+1, k2+1) + ba2(:, k1+l1+1, k2+1));
    elseif l1 == 0           % l2 > 0
      DN_cheb_cos_cos(:,i_lin,j_lin) = 3 * ( ...
        ba2(:, k1+1, abs(k2-l2)+1) + ba2(:, k1+1, k2+l2+1));
    else                     % l1 > 0, l2 > 0  (4 terms)
      DN_cheb_cos_cos(:,i_lin,j_lin) = 3 * ( ...
        ba2(:, abs(k1-l1)+1, abs(k2-l2)+1) + ba2(:, k1+l1+1, abs(k2-l2)+1) ...
      + ba2(:, abs(k1-l1)+1, k2+l2+1)      + ba2(:, k1+l1+1, k2+l2+1));
    end
  end
end

%% Eigenvalues of the linear part (OK-specific)
L1 = intval(L(1)); L2 = intval(L(2));
[kx, ky] = ndgrid(0:m1-1, 0:m2-1);
bkL2 = (kx*L1).^2 + (ky*L2).^2;
mu_k = bkL2 .* (-epsilon^2*bkL2 + 1) - sigma;

%% Midpoint eigendecomposition + Neumann correction
% Linearized operator: A = diag(mu_k) - diag(bkL2) * DN
DN_bar = permute(DN_cheb_cos_cos(1,:,:), [2,3,1]);
A_k0 = diag(mu_k(:)) - bkL2(:) .* DN_bar;
[V, D] = eig(mid(A_k0));
iA_bar = (V * intval(D)) / V;
ilam = intval(diag(D));
iV = intval(V);

iVi=inv(iV);phi=getting_integrated_exponential(real(ilam),h);
integrated_semigroup = abs(iV) * diag(phi) * abs(iVi);

% Correct for time-varying part via Neumann series
% First coefficient: residual from eigendecomposition
DN_cheb_cos_cos(1,:,:) = permute(iA_bar - diag(mu_k(:)) + bkL2(:) .* DN_bar, [3,1,2]);
% Higher coefficients: multiply by bkL2 (derivative factor in A(t) - A_k0)
DN_cheb_cos_cos(2:end,:,:) = permute(bkL2(:), [3,1,2]) .* DN_cheb_cos_cos(2:end,:,:);

% Transform signed defect coefficients before estimating their absolute
% time supremum. The coefficients represent bar A - A(t); its common minus
% sign does not affect the bound for |V^{-1}(A(t)-bar A)|.
delta_eig_sup = intval(zeros(m1m2, m1m2));
for j = 1:m1m2
  defect_cheb = DN_cheb_cos_cos(:,:,j);
  transformed_cheb = (iV \ defect_cheb.').';
  if isreal(transformed_cheb)
    delta_eig_sup(:,j) = chebmag_oneside(transformed_cheb);
  else
    % The extrema routine is for real series. For complex series, retain
    % the cancellation in each transformed coefficient and use |T_n|<=1.
    delta_eig_sup(:,j) = (abs(transformed_cheb(1,:)) ...
      + 2*sum(abs(transformed_cheb(2:end,:)),1)).';
  end
end

integrated_semigroup_alt = abs(iV) * diag(phi);
neumann_defect = intval(max(0,sup(integrated_semigroup_alt * delta_eig_sup)));
if ~(sup(norm(neumann_defect,inf)) < 1)
  error('getting_int_semi_matrix:Neumann', ...
    'Neumann correction could not be validated: ||K||_inf is not below 1.');
end
int_semi_mat = (intval(eye(m1m2)) - neumann_defect) \ integrated_semigroup;

% Version with bkL2 factor on source modes (for Z-bound / nonlinear Lipschitz)
% int_semi_mat_q(k,l) = int_semi_mat(k,l) * bkL2(l)
int_semi_mat_q = int_semi_mat .* bkL2(:).';

if nargout>6
    finite=struct('V',iV,'Vi',iVi,'lambda',ilam,'K',neumann_defect,'h',intval(h), ...
        'Bneu',intval(eye(m1m2))-neumann_defect);
    Sbar=getting_matrix_exp_bound(finite,h);
    finite.error=neumann_defect*(finite.Bneu\Sbar);
end

end
