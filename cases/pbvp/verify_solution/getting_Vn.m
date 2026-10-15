function [Vn,Vnq,e_K_inf,e_inf_K,ba2_norm] = getting_Vn(ba,h,m,params,L,nu)
% params = [lambda,sigma,epsilon];
% lambda  = params(1);
sigma   = params(2);
epsilon = params(3);
%
[~,ba2] = powerconvcos(ba,2);
size_ba2 = size(ba2);
n_ba2 = size_ba2(1); N1_ba2 = size_ba2(2); N2_ba2 = size_ba2(3);
ba2_sup = reshape(chebmag(reshape(ba2,n_ba2,N1_ba2*N2_ba2)),N1_ba2,N2_ba2);
ba2_norm = 3*wnorm(ba2_sup,nu);
%
% setting for OK
% g_ba = 3*ba2_norm; % g(a_bar)
e_K_inf = 3*getting_cm(ba2_sup,m,nu);
e_inf_K = 3*getting_cinf(ba2_sup,m,nu);
ba2_sup(1,1)=0; c_conv=3*wnorm(ba2_sup,nu);
e_K_inf = min(e_K_inf,c_conv);
e_inf_K = min(e_inf_K,c_conv);
%
%% Index k, ell
Chebyshev = size(ba2,1);
m1 = m(1); m2 = m(2); m1m2 = prod(m);
k = (0:m1m2-1)'; ell = 1:m1m2-1;
ba2 = reshape(ba2,Chebyshev,N1_ba2*N2_ba2);
%
%
DN_cheb_cos_cos = intval(zeros(Chebyshev,m1m2,m1m2));
DN_cheb_cos_cos(:,:,1) = 3*ba2(:,k+1);
for j=ell+1
  DN_cheb_cos_cos(:,:,j) = 3*(ba2(:,abs(k-j)+1) + ba2(:,abs(k+j)+1));
end
%
L1 = intval(L(1)); L2 = intval(L(2));
[kx,ky] = ndgrid(0:m1-1,0:m2-1);
% mu_k = -beta1*(k*L).^4 - beta2*(k*L).^2 + 1;
% mu_k = lambda - (1-(kx*L1).^2-(ky*L2).^2).^2 ;
bkL2 = (kx*L1).^2 + (ky*L2).^2;
mu_k = bkL2 .* (-epsilon^2*bkL2 + 1) - sigma;
%
% %% old
% DN_sup = intval(zeros(m1m2,m1m2));
% for j=1:m1m2
%   DN_sup(:,j) = chebrange_oneside(DN_cheb_cos_cos(:,:,j));
% end
% %
% A_sup = diag(mu_k(:)) - bkL2(:) .* DN_sup;
% % ilam = wopnorm(A_sup,nu);
% % Vn = (exp(ilam*h)-1)./(real(ilam));
% [V,D] = eig(mid(A_sup));
% %
% iV = intval(zeros(m1m2));
% ilam = intval(zeros(m1m2,1));
% for j=1:m1m2
%   [ilam(j),iV(:,j)] = verifyeig(A_sup,D(j,j),V(:,j));
% end
% %
% Vn = wopnorm(abs(iV) * diag((exp(ilam*h)-1)./(real(ilam))) * abs(inv(iV)),nu);
% %
%% New
DN_bar = permute(DN_cheb_cos_cos(1,:,:),[2,3,1]);
A_k0 = diag(mu_k(:)) - bkL2(:) .* DN_bar;
[V,D] = eig(mid(A_k0));
iA_bar = (V*intval(D))/V;
ilam = intval(diag(D))
iV = intval(V);
integrated_semigroup = abs(iV) * diag((exp(real(ilam)*h)-1)./(real(ilam))) * abs(inv(iV));
%
DN_cheb_cos_cos(1,:,:) = permute(iA_bar - diag(mu_k(:))+bkL2(:).*DN_bar, [3,1,2]);
DN_cheb_cos_cos(2:end,:,:) = permute(bkL2(:),[3,1,2]) .* DN_cheb_cos_cos(2:end,:,:);
delta_sup = intval(zeros(m1m2,m1m2));
for j=1:m1m2
  delta_sup(:,j) = chebmag_oneside(DN_cheb_cos_cos(:,:,j));
end
integrated_semigroup_alt = abs(iV) * diag((exp(real(ilam)*h)-1)./(real(ilam)));
B = eye(m1m2)-integrated_semigroup_alt*(abs(iV\delta_sup));
Vn  = wopnorm(B\integrated_semigroup,nu);
Vnq = wopnorm((B\integrated_semigroup) .* bkL2(:).',nu);
end
