function [Phi,r1] = verify_variationaleq_KS(ba,h,m,params,L)
% NOTE: For KS, m should be kept small (e.g., m=4-5) because mu_k ~ k^4
% growth makes Z_1 Chebyshev aliasing bound scale as max(|mu_k|)/(2*n_cheb).
% Use more Chebyshev modes (n_pad) if larger m is needed.
% verify_variationaleq_KS  Rigorous verification of variational equation for 1D KS.
%
% Input:
%   ba     - [n, N] intval array of one-sided Chebyshev x sineFourier coefficients
%   h      - time step size
%   m      - truncation size for the variational equation
%   params - [beta1; beta2] physical parameters
%   L      - spatial frequency scaling
%
% Output:
%   Phi - [n, m, m] fundamental matrix solution (Chebyshev coefficients)
%   r1  - rigorous error bound (radii polynomial)

[n, N_total] = size(ba);

% Nop_coeff for KS: -(1/2)(u^2)_x
% Row 1 (j=0): no zero-order derivative terms
% Row 2 (j=1): first derivative, column 3 (i=2, power 2): coefficient -1/2
Nop_coeff = [0, 0, 0; 0, 0, -1/2];

nu_var = 1.0;

n_pad = 30;
ba_ext = zeros(n + n_pad, N_total);
ba_ext(1:n, :) = mid(ba);
ba = ba_ext;
n = n + n_pad;

[Phi, ~, Y0, Z0, Z1] = rig_1D_KS_variational_cheb(ba, m, Nop_coeff, h, params, nu_var, L);

r1 = Y0 ./ (1 - intval(Z0) - intval(Z1));

Phi = reshape(Phi, n, m, m);
end
