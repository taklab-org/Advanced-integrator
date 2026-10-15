function [err,err_at_endpoint,Wh_at_t1,Wh_J,Wm,W_h,ba_X,kappa] = verify_local_existence(eps_all,delta,ba,h,m,params,L,nu)
% params = [lambda,sigma,epsilon];
lambda  = params(1);
sigma   = params(2);
epsilon = params(3);
q = 2; % order of algebraic weights
%
% n = size(ba,1);
[n,N1,N2] = size(ba);
%% estimate the evolution operator
ba_sup = reshape(chebmag(reshape(intval(ba),n,N1*N2)),N1,N2);
ba_X = wnorm(ba_sup,nu,q);
% ba_X_alt = wnorm(permute(sum(abs(intval(ba)),1),[2,3,1]),nu,q);
% ba_Y = wnorm(ba_sup,nu,0);
%
[~,ba2] = powerconvcos(intval(ba),2);
size_ba2 = size(ba2);
n_ba2 = size_ba2(1); N1_ba2 = size_ba2(2); N2_ba2 = size_ba2(3);
ba2_sup = reshape(chebmag(reshape(ba2,n_ba2,N1_ba2*N2_ba2)),N1_ba2,N2_ba2);
% ba2_sup_alt = permute(sum(abs(ba2),1),[2,3,1]);
ba2_X = wnorm(ba2_sup,nu,q);
% ba2_X_alt = wnorm(ba2_sup_alt,nu,q);
ba2_Y = wnorm(ba2_sup,nu,0);
%
bmL2 = (m(1)*L(1)).^2 + (m(2)*L(2)).^2;
mu_m = bmL2 .* (-epsilon^2*bmL2 + 1) - sigma; % eigenvalues
%
% setting for DC
g_ba_X = 3*ba2_X; % g(a_X)
g_ba_Y = 3*ba2_Y; % g(a_Y)
Lq_hat = sup(max(intval(L))^q);
cmq = Lq_hat*getting_cm(ba2_sup,m,nu,q);
cm0 = Lq_hat*getting_cm(ba2_sup,m,nu,0);
cinfq = Lq_hat*getting_cinf(ba2_sup,m,nu,q);
%
ba2_sup(1,1)=0; c_conv=3*Lq_hat*wnorm(ba2_sup,nu,q);
cmq = min(cmq,c_conv);
cm0 = min(cm0,3*Lq_hat*wnorm(ba2_sup,nu,0));
cinfq = min(cinfq,c_conv);
% cinf0 = Lq_hat*getting_cinf(ba2_sup,m,nu,0);
%
gamma = (1 + m(1)*L(1) + m(2)*L(2)).^q; % <m>^q
vtheta = mu_m+gamma*Lq_hat*g_ba_X;
vtheta0 = mu_m+gamma*Lq_hat*g_ba_Y;
if vtheta>100
  disp('The alpha value is too large. Take more ''m''.')
  Wh_at_t1 = NaN; W_h = NaN; Wh_J = NaN; err = NaN; err_at_endpoint=NaN; Wm=NaN; kappa=NaN;
  return
end
%
% Bounds for the evolution operator in the tail
W_infq = gamma*(exp(vtheta*h)-1)/vtheta;
W_inf0 = (exp(vtheta0*h)-1)/vtheta0;
W_bar_infq = (W_infq-gamma*h)/vtheta;
W_bar_inf0 = (W_inf0-h)/vtheta0;
W_inf_supq = gamma*max(1,exp(vtheta*h));
W_inf_sup0 = max(1,exp(vtheta0*h));

%% solve variational problem
Nop_coeff = [0,0,0,0;0,0,0,0;0,0,0,-1];
nu_var = 1.4;
m1m2 = prod(m);
%
% padding chebyshev for variational problems
n_pad = 15;
ba_ext = zeros(n+n_pad,N1,N2);
ba_ext(1:n,:,:) = ba;
ba = ba_ext;
n = n + n_pad;
%
[Phi,~,Y0,Z0,Z1] = rig_2D_diblock_copolymer_variational(ba,m,Nop_coeff,h,epsilon,sigma,nu_var,L,lambda);
r1 = Y0./(1-intval(Z0)-intval(Z1));
Phi = reshape(Phi,n,m1m2,m1m2);
%
nu_var = 1.65;
[Psi,~,Y0,Z0,Z1] = rig_2D_diblock_copolymer_adjoint(ba,m,Nop_coeff,h,epsilon,sigma,nu_var,L,lambda);
r2 = Y0./(1-intval(Z0)-intval(Z1));
Psi = reshape(Psi,n,m1m2,m1m2);
%
% A technique to estimate a uniform bound using interval arithmetic
[Wm,Wm_J,Wm_at_t1,Wm_bar] = getting_Wm(Phi,Psi,r1,r2,h,n,m,L,nu,q);
Wm = sup(Wm); Wm_J = sup(Wm_J);
Phi(2:end,:,:) = 2*Phi(2:end,:,:); % convert to two-sided chebyshev
Wm_at_t1 = min(wopnorm(reshape(sum(intval(Phi),1), m1m2, m1m2)+r1, nu), Wm_at_t1);
%
kappa = 1 - Wm*W_bar_infq*cmq*cinfq;
kappa = inf(kappa);
%
if kappa>0
  U_matrix = [Wm, Wm*W_infq*cmq;...
    Wm*W_infq*cinfq, W_inf_supq]/kappa;
  U_matrix_at_t1 = [Wm_at_t1+Wm_J*W_bar_infq*cmq*cinfq*U_matrix(1,1),...
    Wm_J*cmq*(W_infq+W_bar_infq*cinfq*U_matrix(1,2));...
    Wm*cinfq*(W_inf0+W_bar_inf0*cmq*U_matrix(2,1)),...
    exp(vtheta0*h)+Wm*W_bar_inf0*cmq*cinfq*U_matrix(2,2)];
  U_matrix_at_t1_alt = [Wm_at_t1+Wm_bar*W_bar_inf0*cm0*cinfq*U_matrix(1,1),...
    Wm_bar*cm0*(W_inf0+W_bar_inf0*cinfq*U_matrix(1,2));...
    Wm*cinfq*(W_inf0+W_bar_inf0*cmq*U_matrix(2,1)),...
    exp(vtheta0*h)+Wm*W_bar_inf0*cmq*cinfq*U_matrix(2,2)];
  U_matrix_J = [Wm_J*(1+W_bar_infq*cmq*cinfq*U_matrix(1,1)),...
    Wm_J*cmq*(W_infq+W_bar_infq*cinfq*U_matrix(1,2));...
    Wm*cinfq*(W_inf0+W_bar_inf0*cmq*U_matrix(2,1)),...
    W_inf_sup0+Wm*W_bar_inf0*cmq*cinfq*U_matrix(2,2)];
  U_matrix_J_alt = [Wm_J+Wm_bar*W_bar_inf0*cm0*cinfq*U_matrix(1,1),...
    Wm_bar*cm0*(W_inf0+W_bar_inf0*cinfq*U_matrix(1,2));...
    Wm*cinfq*(W_inf0+W_bar_inf0*cmq*U_matrix(2,1)),...
    W_inf_sup0+Wm*W_bar_inf0*cmq*cinfq*U_matrix(2,2)];
  W_h = sup(norm(U_matrix,1)); % Matrix 1-norm
  ba_X = sup(ba_X);
  L_rho = @(x) 3 * (2*ba_X+x) .* x; % Special for cubic nolinearity
  F = @(x) W_h*(eps_all+h*(Lq_hat * L_rho(x) .* x + delta))-x;
  Wh_at_t1 = min([W_h,sup(norm(U_matrix_at_t1,1)),sup(norm(U_matrix_at_t1_alt,1))]);
  Wh_J = min([W_h,sup(norm(U_matrix_J,1)),sup(norm(U_matrix_J_alt,1))]);
else
  disp('Linearized problem is not solved (kappa<0)')
  err_at_endpoint=NaN; Wm=NaN;
  Wh_at_t1 = NaN; W_h = NaN; Wh_J = NaN; err = NaN;
  return
end
%
%
%% verify the contraction mapping:
[xx,xx_cand,data] = verifynlssall(F,infsup(0,1));
if ~any(xx>0)
  disp('contraction mapping is not verified!')
  err = NaN;
  return
end
while(1)
  if isempty(xx_cand)
    err = sup(xx(:,all(F((1+eps)*sup(xx))<0,1)));
    break
  else
    [xx,xx_cand,data] = verifynlssall(data);
  end
end
%% error bound at the endpoint
err_at_endpoint = intval(Wh_at_t1)*eps_all+intval(Wh_J)*h*(Lq_hat*L_rho(err)*err+intval(delta));
err_at_endpoint = min(err,sup(err_at_endpoint));
disp(['error at endpoint = ',num2str(err_at_endpoint)])
%
