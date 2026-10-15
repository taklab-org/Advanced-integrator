function [Phi,r1,verification] = verify_variationaleq(ba,h,m,params,L,n_pad_in)
% params = [lambda,sigma,epsilon];
lambda  = params(1);
sigma   = params(2);
epsilon = params(3);
%
[n_orig,N1,N2] = size(ba);
assert(all(rad(intval(ba))==0,'all'),'GE:VariationalCenter', ...
  'Variational input must be the exact stored polynomial, not an uncertain coefficient box.');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve variational problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nop_coeff = [0,0,0,0;0,0,0,0;0,0,0,-1];
nu_var = 1.0;
m1m2 = prod(m);
%
% padding chebyshev for variational problems
if exist('n_pad_in','var') && ~isempty(n_pad_in)
  n_pad = n_pad_in;
else
  n_pad = 30;  % m=[7,7] adequate
end
ba_ext = zeros(n_orig+n_pad, N1, N2);
ba_ext(1:n_orig,:,:) = mid(ba);
ba = ba_ext;
n = n_orig + n_pad;
%
[Phi,~,Y0,Z0,Z1] = rig_2D_diblock_copolymer_variational_cheb(ba,m,Nop_coeff,h,epsilon,sigma,nu_var,L,lambda);
den = 1-intval(Z0)-intval(Z1);
assert(all(inf(den)>0,'all'),'GE:VariationalContraction','Z0+Z1 is not strictly below one.');
r1 = intval(Y0)./den;
assert(all(isfinite(sup(r1)) & inf(r1)>=0,'all'),'GE:VariationalRadius','Invalid variational enclosure radius.');
verification=struct('passed',true,'Y0',intval(Y0),'Z0',intval(Z0), ...
  'Z1',intval(Z1),'margin',den,'n_cheb',n,'m',m,'h',h);
%
Phi = reshape(Phi,n,m1m2,m1m2);
