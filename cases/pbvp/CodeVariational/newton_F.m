function [c,k] = newton_F(mu_k,Nop_coeff,a_bar,c,h,e_j,L,lambda)
% Newton iteration using fused DF_lin_fin_dim (block + cache).
if ~exist('lambda','var')
    lambda = [];
end
tol = 5e-11;
size_c = size(c);
F = F_lin_fin_dim(mu_k,Nop_coeff,a_bar,c,h,e_j,L,lambda);
F = reshape(F,[numel(c),1]);
k_max = 100;
k=0;
while k<=k_max && norm(F,1)> tol
    DF = DF_lin_fin_dim(mu_k,Nop_coeff,a_bar,h,c,L,lambda);
    DF = reshape(DF,[numel(c),numel(c)]);
    c = reshape(c,[numel(c),1]);
    c = c - DF\F;
    c = reshape(c,size_c );
    F = F_lin_fin_dim(mu_k,Nop_coeff,a_bar,c,h,e_j,L,lambda);
    F = reshape(F,[numel(c),1]);
    k=k+1;
end
if norm(F)> tol
    error('Did not converged')
end
if isnan(norm(F,1))
    error('Did not converged')
end
display(['||F(a)|| = ',num2str(norm(F,1)),', # Newton iterations = ',num2str(k)])

end
