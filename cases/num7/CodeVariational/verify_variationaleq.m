function [Phi,r1] = verify_variationaleq(ba,h,m,params,L)
% verify_variationaleq  Rigorous verification of the variational equation
%   [Phi,r1] = verify_variationaleq(ba,h,m,params,L)
%
% Input:
%   ba     - [n, N1, N2] intval array of Chebyshev-Fourier coefficients
%   h      - time step size
%   m      - [m1, m2] truncation sizes for the variational equation
%   params - [lambda, sigma, epsilon] physical parameters
%   L      - [L1, L2] spatial frequency scaling
%
% Output:
%   Phi - [n, m1*m2, m1*m2] fundamental matrix solution
%   r1  - rigorous error bound (radii polynomial)

lambda  = params(1);
sigma   = params(2);
epsilon = params(3);

[n,N1,N2] = size(ba);

Nop_coeff = [0,0,0,0;0,0,0,0;0,0,0,-1];

nu_var = 1.2;

m1m2 = prod(m);

n_pad = 28;
ba_ext = zeros(n+n_pad,N1,N2);
ba_ext(1:n,:,:) = mid(ba);
ba = ba_ext;
n = n + n_pad;

[Phi,~,Y0,Z0,Z1] = rig_2D_diblock_copolymer_variational(ba,m,Nop_coeff,h,epsilon,sigma,nu_var,L,lambda);

r1 = Y0./(1-intval(Z0)-intval(Z1));

Phi = reshape(Phi,n,m1m2,m1m2);
end
