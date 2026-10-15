function [F] = F_lin_fin_dim_KS(mu_k,Nop_coeff,a_bar,c,h,e_j,L,lambda)
% F_lin_fin_dim_KS  Variational equation zero map for KS (sine basis).
%   Same as F_lin_fin_dim but calls N_op_lin_KS (convapbqsin).

if ~exist('lambda','var')
    lambda = [];
end

Mc = size(c);
Ma = size(a_bar);
c_temp = c;
mu_k_temp = mu_k;
c = zeros(Ma);
mu_k = zeros(Ma);
for i = 1:length(Ma)
    otherdims_a{i} = 1:Ma(i);
    otherdims_c{i} = 1:Mc(i);
    otherdims_zeros{i} = 1;
end
c(otherdims_c{:}) = c_temp;
mu_k(otherdims_c{:}) = mu_k_temp;
N = size(a_bar)-1;

a_bar_pad = zeros(Ma+1);
a_bar_pad(otherdims_a{:}) = a_bar;

c_pad = zeros(Ma+1);
c_pad(otherdims_a{:}) = c;

mu_k_pad = zeros(Ma+1);
mu_k_pad(otherdims_a{:}) = mu_k;

nop = N_op_lin_KS(Nop_coeff,a_bar_pad,c_pad);
M = size(nop)-1;
dim = length(M)-1;
s = arrayfun(@(k) 0:k, M(1:end-1), 'UniformOutput', false);
k = cell(1, dim);
[k{:}] = ndgrid(s{:});
if (size(Nop_coeff,1)-1) > 0
    K_sq = zeros(size(k{1}));
    for ell = 1:dim
        K_sq = K_sq + (k{ell}*L(ell)).^(size(Nop_coeff,1)-1);
    end
else
    K_sq = ones(size(k{1}));
end
N_rep = ones(size(M));
N_rep(end) = M(end)+1;
K_sq = repmat(K_sq,N_rep);
nop = K_sq.*nop;

psi = h/2*(mu_k_pad.*c_pad + nop);
F = lambda_op(c_pad) + T_op(psi);
F = F(otherdims_c{:});

F(otherdims_c{1:(end-1)}, 1) = -e_j;
for j = 0:N(end)
    if j == 0
        F(otherdims_c{1:end-1}, 1) = F(otherdims_c{1:end-1}, 1) + c(otherdims_c{1:end-1}, j+1);
    else
        F(otherdims_c{1:end-1}, 1) = F(otherdims_c{1:end-1}, 1) + 2*((-1)^j)*c(otherdims_c{1:end-1}, j+1);
    end
    if ~isempty(lambda) && j > 0
        F(otherdims_zeros{1:end-1},j+1) = 2*j*c(otherdims_zeros{1:end-1},j+1);
    end
end

end
