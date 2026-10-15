function [norm_B] = norm_B_nu_int(B, nu, N)
% Vectorized matrix norm (from V2's norm_B_nu_int_v2).
% Replaces per-row loop with single matrix-vector multiply.

dims = N + 1;
dim = length(N);
numel_c = prod(dims);

B = reshape(abs(B), [numel_c, numel_c]);

% Build nu_vect using same tensor product approach as norm_1_nu_int
w1 = 2 * intval(ones(dims(1), 1));
w1(1) = intval(1);
nu_vect = w1;

for d = 2:dim-1
    w = 2 * intval(ones(1, dims(d)));
    w(1) = intval(1);
    shape = ones(1, dim);
    shape(d) = dims(d);
    nu_vect = nu_vect .* reshape(w, shape);
end

if dim >= 2
    w_cheb = intval(zeros(1, dims(end)));
    for k = 0:dims(end)-1
        if k == 0
            w_cheb(k+1) = intval(1);
        else
            w_cheb(k+1) = 2 * nu^k;
        end
    end
    shape = ones(1, dim);
    shape(dim) = dims(dim);
    nu_vect = nu_vect .* reshape(w_cheb, shape);
end

nu_vect_flat = reshape(nu_vect, [numel_c, 1]);

% Single matrix-vector multiply replaces numel_c calls to norm_1_nu_int
vect = B * nu_vect_flat;

% Inverse nu weights for the max operation
nu_vect_inv = intval(zeros(dims));
otherdims = repmat({':'}, 1, dim-1);
for i = 0:N(end)
    nu_vect_inv(otherdims{:}, i+1) = nu^(-i);
end

norm_B = reshape(vect, dims) .* nu_vect_inv;

for i = 1:dim
    norm_B = max(norm_B);
end

end
