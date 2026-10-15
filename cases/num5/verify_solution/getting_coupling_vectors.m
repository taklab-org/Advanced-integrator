function [psi_vec, c_out_per_mode] = getting_coupling_vectors(ba2, m, ~, nu, coefficient_mode)
% Mode-resolved convolution coupling in ndgrid column-major ordering.
% The third argument is retained for existing callers; the support is taken
% from the full coefficient array ba2, not from the original truncation.
%
% Default 'majorant': ba2 contains entrywise time-supremum bounds for |u^2|.
% 'signed': ba2 encloses the signed coefficients of a time-independent u^2.
% In signed mode, reflected coefficients are added before taking abs so that
% cancellation is retained. Neither output includes the factor 3 or bkL2.

if nargin < 5 || isempty(coefficient_mode)
  coefficient_mode = 'majorant';
end
switch coefficient_mode
  case 'majorant'
    kernel = abs(intval(ba2));
  case 'signed'
    kernel = intval(ba2);
  otherwise
    error('getting_coupling_vectors:Mode', ...
      'coefficient_mode must be ''majorant'' or ''signed''.');
end

m1 = m(1); m2 = m(2); m1m2 = prod(m);
[B1, B2] = size(kernel);
[kx_grid, ky_grid] = ndgrid(0:m1-1, 0:m2-1);
kx_fin = kx_grid(:); ky_fin = ky_grid(:);

% A sentinel row and column represent exact zero outside the kernel support.
kernel_ext = intval(zeros(B1+1, B2+1));
kernel_ext(1:B1, 1:B2) = kernel;

%% Tail -> finite: maximum weighted column contribution to each finite row.
psi_vec = intval(zeros(m1m2, 1));
for k_lin = 1:m1m2
  k1 = kx_fin(k_lin); k2 = ky_fin(k_lin);
  % A nonzero coefficient requires ell_j <= (B_j-1) + k_j.
  [ell1_grid, ell2_grid] = ndgrid(0:B1+k1-1, 0:B2+k2-1);
  in_Fm = (ell1_grid < m1) & (ell2_grid < m2);
  ell1 = ell1_grid(~in_Fm); ell2 = ell2_grid(~in_Fm);
  if isempty(ell1), continue; end

  coefficient = reflected_coefficient(kernel_ext, k1, k2, ell1, ell2);
  coefficient_abs = abs(coefficient);
  w_ell = mode_weight(ell1, ell2, nu);
  psi_vec(k_lin) = max(mag(coefficient_abs ./ w_ell));
end

%% Finite -> tail: sum all weighted output modes for each finite column.
c_out_per_mode = intval(zeros(m1m2, 1));
for l_lin = 1:m1m2
  ell1 = kx_fin(l_lin); ell2 = ky_fin(l_lin);
  % Multiplication by mode ell extends the output support beyond that of u^2.
  [k1_grid, k2_grid] = ndgrid(0:B1+ell1-1, 0:B2+ell2-1);
  in_Fm = (k1_grid < m1) & (k2_grid < m2);
  k1 = k1_grid(~in_Fm); k2 = k2_grid(~in_Fm);
  if isempty(k1), continue; end

  coefficient = reflected_coefficient(kernel_ext, k1, k2, ell1, ell2);
  coefficient_abs = abs(coefficient);
  w_k = mode_weight(k1, k2, nu);
  c_out_per_mode(l_lin) = sum(coefficient_abs .* w_k);
end
end

function coefficient = reflected_coefficient(kernel, k1, k2, ell1, ell2)
% Sum over the distinct signed source modes. An axis with ell_j=0 has only
% one source sign; repeated coefficient indices for k_j=0 still contribute.
coefficient = kernel_entry(kernel, abs(k1-ell1), abs(k2-ell2));
coefficient = coefficient + double(ell1>0) .* ...
  kernel_entry(kernel, k1+ell1, abs(k2-ell2));
coefficient = coefficient + double(ell2>0) .* ...
  kernel_entry(kernel, abs(k1-ell1), k2+ell2);
coefficient = coefficient + double((ell1>0) & (ell2>0)) .* ...
  kernel_entry(kernel, k1+ell1, k2+ell2);
end

function value = kernel_entry(kernel, j1, j2)
% Clamp only to the zero sentinel, never to a coefficient at the boundary.
i1 = min(j1+1, size(kernel,1));
i2 = min(j2+1, size(kernel,2));
value = kernel(sub2ind(size(kernel), i1, i2));
end

function w = mode_weight(k1, k2, nu)
alpha = 2.^(double(k1>0) + double(k2>0));
w = alpha .* intval(nu).^(k1+k2);
end
