function [F] = F_lin_fin_dim_int(mu_k,Nop_coeff,a_bar,c,h,e_j,L,lambda,a_bar_pow,nop_scalar,support_sz)
% Optional args 9-11: precomputed a_bar^(P-1), scalar coefficient, support size
% for direct convolution (Improvements A+B). If not provided, falls back to FFT.
if ~exist('lambda','var')
    lambda = [];
end
if ~exist('a_bar_pow','var')
    a_bar_pow = [];
end
Mc = size(c);
Ma = size(a_bar);
c_temp = c;
mu_k_temp = mu_k;
c = intval(zeros(Ma));
mu_k = intval(zeros(Ma));
for i = 1:length(Ma)
    otherdims_a{i} = 1:Ma(i);
    otherdims_c{i} = 1:Mc(i);
    otherdims_zeros{i} = 1;
end
c(otherdims_c{:}) = c_temp;
mu_k(otherdims_c{:}) = mu_k_temp;
N = size(a_bar)-1;

a_bar_pad = intval(zeros(Ma+1));
a_bar_pad(otherdims_a{:}) = a_bar;

c_pad = intval(zeros(Ma+1));
c_pad(otherdims_a{:}) = c;

mu_k_pad = intval(zeros(Ma+1));
mu_k_pad(otherdims_a{:}) = mu_k;

if ~isempty(a_bar_pow)
    % [Opt4] Fast path: direct convolution using precomputed a_bar^(P-1)
    nop = nop_scalar * cosconv_small_support(a_bar_pow, c_pad, support_sz, size(a_bar_pad));
else
    nop = N_op_lin_int(Nop_coeff,a_bar_pad,c_pad);
end
M = size(nop)-1;
dim = length(M)-1;
s = arrayfun(@(k) 0:k, M(1:end-1), 'UniformOutput', false);
k = cell(1, dim);
[k{:}] = ndgrid(s{:});
if (size(Nop_coeff,1)-1) > 0
    K_sq = intval(zeros(size(k{1})));
    for ell = 1:dim
        K_sq = K_sq + (k{ell}*L(ell)).^(size(Nop_coeff,1)-1);
    end
else
    K_sq = intval(ones(size(k{1})));
end
N_rep = ones(size(M));
N_rep(end) = M(end)+1;
K_sq = repmat(K_sq,N_rep);
nop = K_sq.*nop;

psi = h/2*(mu_k_pad.*c_pad + nop);
F = lambda_op_int(c_pad) + T_op_int(psi);
F = F(otherdims_c{:});

F(otherdims_c{1:end-1}, 1) = F(otherdims_c{1:end-1}, 1) - e_j + c(otherdims_c{1:end-1}, 1);
for j = 1:N(end)
    F(otherdims_c{1:end-1}, 1) = F(otherdims_c{1:end-1}, 1) + 2*((-1)^j)*c(otherdims_c{1:end-1}, j+1);
end
if ~isempty(lambda)
    for j = 1:N(end)
        F(otherdims_zeros{1:end-1},j+1) = 2*j*c(otherdims_zeros{1:end-1},j+1);
    end
end
end
