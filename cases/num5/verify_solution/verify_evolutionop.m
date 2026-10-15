function [W_tau, W_tau_q, W_J, W_J_q, W_at_t1, ba_norm, kappa, gam] = verify_evolutionop(ba,tau,m,params,L,nu,mu_ast)
% params = [lambda,sigma,epsilon];
lambda  = params(1);
sigma   = params(2);
epsilon = params(3);
% q = 2; % order of algebraic weights
%
% n = size(ba,1);
[n,N1,N2] = size(ba);
%% estimate the evolution operator
ba_norm = wnorm(reshape(chebmag(reshape(intval(ba),n,N1*N2)),N1,N2),nu);
%
[~,ba2] = powerconvcos(intval(ba),2);
size_ba2 = size(ba2);
n_ba2 = size_ba2(1); N1_ba2 = size_ba2(2); N2_ba2 = size_ba2(3);
ba2_sup = reshape(chebmag(reshape(ba2,n_ba2,N1_ba2*N2_ba2)),N1_ba2,N2_ba2);
% ba2_sup_alt = permute(sum(abs(ba2),1),[2,3,1]);
ba2_norm = wnorm(ba2_sup,nu);
%
%
% setting for DC
g_ba = 3*ba2_norm; % g(a_X)
e_1_inf = getting_cm(ba2_sup,m,nu);
e_inf_1 = getting_cinf(ba2_sup,m,nu);
%
ba2_sup(1,1)=0; c_conv=3*wnorm(ba2_sup,nu);
e_1_inf = min(e_1_inf,c_conv);
e_inf_1 = min(e_inf_1,c_conv);
%
xi = intval(.8); % optimal?
gam = intval(0.5); % should be greater than 1/2
Cinf = compute_Cinf(params,m,L,gam,xi);
iota = (1-xi)*mu_ast;
beta_one_minus_gam = beta(1-gam, 1-gam);
vtheta = Cinf * g_ba * beta_one_minus_gam;
tvtheta = Cinf * g_ba / (1-gam);
%
% Bounds for the evolution operator : W_tau_q
W_infinite_sup_q = Cinf * exp(vtheta*tau^(1-gam));
W_infinity_q = tau / (1-gam) * W_infinite_sup_q;
W_bar_infty_q = W_infinity_q/2 * tau^(1-gam) * beta_one_minus_gam;
%
% Bounds for the evolution operator : W_tau
W_infinite_sup = exp(tvtheta * tau^(1-gam));
if iota==0
  W_infinity = tau * W_infinite_sup;
else
  W_infinity = (1 - exp(-iota * tau))/iota * W_infinite_sup;
end
W_hat_infty = tau^(1-gam) / (1-gam) * W_infinite_sup_q;
W_bar_infty = tau^(2-gam) / (1-gam) / (2-gam) * W_infinite_sup_q;
% ;
% return
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
[Phi,~,Y0,Z0,Z1] = rig_2D_diblock_copolymer_variational(ba,m,Nop_coeff,tau,epsilon,sigma,nu_var,L,lambda);
r1 = Y0./(1-intval(Z0)-intval(Z1));
Phi = reshape(Phi,n,m1m2,m1m2);
%
nu_var = 1.5;
[Psi,~,Y0,Z0,Z1] = rig_2D_diblock_copolymer_adjoint(ba,m,Nop_coeff,tau,epsilon,sigma,nu_var,L,lambda);
r2 = Y0./(1-intval(Z0)-intval(Z1));
Psi = reshape(Psi,n,m1m2,m1m2);
%
% A technique to estimate a uniform bound using interval arithmetic
[Wmq,Wm0,Wmq_J,Wm0_J,Wm_at_t1] = getting_Wm(Phi,Psi,r1,r2,tau,n,m,L,nu);
Wmq = sup(Wmq); Wm0 = sup(Wm0); Wmq_J = sup(Wmq_J);
Phi(2:end,:,:) = 2*Phi(2:end,:,:); % convert to two-sided chebyshev
Wm_at_t1 = min(wopnorm(reshape(sum(intval(Phi),1), m1m2, m1m2)+r1, nu), Wm_at_t1);
%
kappa = 1 - Wmq*W_bar_infty_q*e_1_inf*e_inf_1;
tkappa = 1 - Wmq*W_bar_infty*e_1_inf*e_inf_1;
kappa = inf(kappa);
%
if kappa>0
  U_matrix = [tau^gam * Wmq, Wmq*W_infinity_q*e_1_inf;...
    Wmq*W_infinity_q*e_inf_1, W_infinite_sup_q]/kappa;
  W_tau_q = (norm(U_matrix,1)); % Matrix 1-norm
  %
  tU_matrix = [Wm0, Wmq*W_infinity*e_1_inf;...
    Wmq*W_hat_infty*e_inf_1, W_infinite_sup]/tkappa;
  W_tau = (norm(tU_matrix,1)); % Matrix 1-norm
  %
  U_matrix_at_t1 = [Wm_at_t1+Wmq_J*W_bar_infty*e_1_inf*e_inf_1*tU_matrix(1,1),...
    Wmq_J*e_1_inf*(W_infinity+W_bar_infty*e_inf_1*tU_matrix(1,2));...
    e_inf_1*(Wm0*W_hat_infty+Wmq*W_bar_infty*e_1_inf*tU_matrix(2,1)),...
    exp(-(1-xi)*mu_ast*tau + Cinf*g_ba*tau^(1-gam)/(1-gam))+Wmq*W_bar_infty*e_1_inf*e_inf_1*tU_matrix(2,2)];
  W_at_t1 = norm(U_matrix_at_t1,1);
  %
  U_matrix_J_q = [Wmq_J*(tau^gam + W_bar_infty_q*e_1_inf*e_inf_1*U_matrix(1,1)),...
    Wmq_J*e_1_inf*(W_infinity_q+W_bar_infty_q*e_inf_1*U_matrix(1,2));...
    Wmq*e_inf_1*(W_infinity_q+W_bar_infty_q*e_1_inf*U_matrix(2,1)),...
    W_infinite_sup_q+Wmq*W_bar_infty_q*e_1_inf*e_inf_1*U_matrix(2,2)];
  W_J_q = norm(U_matrix_J_q,1);
  %
  U_matrix_J = [Wm0_J + Wmq_J*W_bar_infty*e_1_inf*e_inf_1*tU_matrix(1,1),...
    Wmq_J*e_1_inf*(W_infinity+W_bar_infty*e_inf_1*tU_matrix(1,2));...
    e_inf_1*(W_hat_infty*Wm0 + Wmq*W_bar_infty*e_1_inf*tU_matrix(2,1)),...
    W_infinite_sup+Wmq*W_bar_infty*e_1_inf*e_inf_1*tU_matrix(2,2)];
  W_J = norm(U_matrix_J,1);
  %
  ba_norm = sup(ba_norm);% ba2_X = sup(ba2_X);
  W_at_t1 = min(W_tau,W_at_t1);
  W_J_q = min(W_tau_q,W_J_q);
  W_J = min(W_tau,W_J);
else
  disp('Linearized problem is not solved (kappa<0)')
  W_tau = NaN; W_tau_q = NaN; W_J = NaN; W_J_q = NaN; W_at_t1 = NaN;
  return
end
