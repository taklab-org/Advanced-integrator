function [Psi,r,Y_0,Z_0,Z_1] = rig_1D_SH_variational(ba,m,Nop_coeff,h,params,nu,L)
% Fused: V2 block interval DF + New1 bounds precomputation + V2 vectorized norms.
a_bar = ba.';
c = zeros([m,size(ba,1)]);
mu_k = mu_k_1DSH(size(c)-1,params,L);
Y_0 = zeros(1,m);
Psi = zeros(numel(c),m);
DF_nm = reshape(DF_lin_fin_dim_int(intval(mu_k),intval(Nop_coeff),intval(a_bar),intval(h),intval(c),intval(L)),[numel(c),numel(c)]);
bounds_ctx = prepare_bounds_context(intval(mu_k),intval(Nop_coeff),intval(a_bar),intval(h),intval(nu),DF_nm,intval(L));
for i = 1:m
  e_j = zeros(m,1);
  e_j(i) = 1;
  [c,~] = newton_F(mu_k,Nop_coeff,a_bar,c,h,e_j,L);
  if i ==1
    [r, Y_0(1,i),Z_0,Z_1] = eval_full_bounds(bounds_ctx,intval(c),intval(e_j));
  else
    Y_0(1,i) = eval_bound_y0(bounds_ctx,intval(c),intval(e_j));
  end
  Psi(:,i) = reshape(c.',[],1);
end
end
