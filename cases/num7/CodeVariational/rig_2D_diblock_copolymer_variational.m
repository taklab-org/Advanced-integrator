function [Psi,r,Y_0,Z_0,Z_1] = rig_2D_diblock_copolymer_variational(ba,m,Nop_coeff,h,epsilon,sigma,nu,L,lambda)
% Fused: V2 block interval DF + New1 bounds precomputation + V2 vectorized norms.
if ~exist('lambda','var')
    lambda = [];
end
a_bar = permute(ba,[2,3,1]);
c = zeros([m,size(ba,1)]);
mu_k = mu_k_dcCH(size(c)-1,epsilon,sigma,L);
Y_0 = zeros(1,m(1)*m(2));
Psi = zeros(numel(c),m(1)*m(2));
DF_nm = reshape(DF_lin_fin_dim_int(intval(mu_k),intval(Nop_coeff),intval(a_bar),intval(h),intval(c),intval(L),lambda),[numel(c),numel(c)]);
bounds_ctx = prepare_bounds_context(intval(mu_k),intval(Nop_coeff),intval(a_bar),intval(h),intval(nu),DF_nm,intval(L),lambda);
for i = 1:m(1)*m(2)
        e_j = zeros(m);
        e_j(i) = 1;
        [c,~] = newton_F(mu_k,Nop_coeff,a_bar,c,h,e_j,L,lambda);
        if i ==1
            [r, Y_0(1,i),Z_0,Z_1] = eval_full_bounds(bounds_ctx,intval(c),intval(e_j));
        else
            Y_0(1,i) = eval_bound_y0(bounds_ctx,intval(c),intval(e_j));
        end
        Psi(:,i) = reshape(permute(c,[3,1,2]),[],1);
end
end
