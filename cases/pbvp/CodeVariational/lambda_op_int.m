function [Lambda] = lambda_op_int(a)
% [Opt6] Vectorized: replace loop with repmat + element-wise multiply.
N = size(a)-1;
weights = 2*(0:N(end));
w_shape = ones(1, length(N));
w_shape(end) = N(end)+1;
w_rep = N+1;
w_rep(end) = 1;
Lambda = repmat(reshape(intval(weights), w_shape), w_rep) .* a;
end
