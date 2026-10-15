function [Phi,r1] = verify_variationaleq(ba,h,m,params,L)
% params = [lambda,sigma,epsilon];
lambda  = params(1);
sigma   = params(2);
epsilon = params(3);
%
[n,N1,N2] = size(ba);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve variational problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nop_coeff = [0,0,0,-1];
Nop_coeff = [0,0,0,0;0,0,0,0;0,0,0,-1];
% nu_var = 1.0;
nu_var = 1.2;
% n = size(ba,1);
m1m2 = prod(m);
%
% padding chebyshev for variational problems
% n_pad = 0;
n_pad = 20;
ba_ext = zeros(n+n_pad,N1,N2);
ba_ext(1:n,:,:) = mid(ba);
ba = ba_ext;
n = n + n_pad;
%
% [Phi,~,Y0,Z0,Z1] = rig_2D_SH_variational(ba,m,Nop_coeff,h,lambda,nu_var,L);
% r1 = Y0./(1-intval(Z0)-intval(Z1));
[Phi,~,Y0,Z0,Z1] = rig_2D_diblock_copolymer_variational(ba,m,Nop_coeff,h,epsilon,sigma,nu_var,L,lambda);
r1 = Y0./(1-intval(Z0)-intval(Z1));
%
Phi = reshape(Phi,n,m1m2,m1m2);
