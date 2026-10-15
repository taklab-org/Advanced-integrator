function [mu_k] = mu_k_dcCH(N,epsilon,sigma,L)
dim = length(N)-1;
s = arrayfun(@(k) 0:k, N(1:end-1), 'UniformOutput', false);
k = cell(1, dim);
[k{:}] = ndgrid(s{:});
K_sq = zeros(size(k{1}));
for i = 1:dim
    K_sq = K_sq + (k{i}*L(i)).^2;
end
mu_k = (K_sq).*(-epsilon^2*K_sq +1)-sigma;
N_rep = ones(size(N));
N_rep(end) = N(end)+1;
mu_k = repmat(mu_k,N_rep);
end
