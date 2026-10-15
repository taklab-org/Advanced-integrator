function [mu_k] = mu_k_1DKS(N,params,L)
% mu_k_1DKS  Linear decay rates for 1D KS in sine basis.
%   mu_k = -beta1*(k*L)^4 - beta2*(k*L)^2  (no +1 term, unlike SH)
%   Uses 0-indexed convention: k = 0:N(1), with mu_0 = 0.
beta1 = params(1);
beta2 = params(2);
dim = length(N)-1;
s = arrayfun(@(k) 0:k, N(1:end-1), 'UniformOutput', false);
k = cell(1, dim);
[k{:}] = ndgrid(s{:});
mu_k = -beta1*(k{1}*L).^4 - beta2*(k{1}*L).^2;
N_rep = ones(size(N));
N_rep(end) = N(end)+1;
mu_k = repmat(mu_k,N_rep);
end
