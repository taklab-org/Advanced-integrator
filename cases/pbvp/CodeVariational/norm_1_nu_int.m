function [n] = norm_1_nu_int(a, nu, nu_vect_precomp)
% Vectorized weighted 1-norm (from V2's norm_1_nu_int_v2).
% Computes weighted 1-norm: sum(|a| .* nu_vect) where
%   nu_vect(k1,...,kd) = prod(2.^(ki>0)) * nu^(k_last)
%
% Optional 3rd argument: precomputed nu_vect (for reuse in norm_B_nu_int)

N = size(a) - 1;
dims = N + 1;
dim = length(N);

if nargin >= 3 && ~isempty(nu_vect_precomp)
    nu_vect = nu_vect_precomp;
else
    % Build nu_vect via tensor product of 1D weight vectors
    % Each spatial dimension d: w(0)=1, w(k)=2 for k>0
    % Last dimension (Chebyshev): w(0)=1, w(k)=2*nu^k for k>0

    % Start with the first dimension
    w1 = 2 * intval(ones(dims(1), 1));
    w1(1) = intval(1);
    nu_vect = w1;

    for d = 2:dim-1
        % Spatial dimensions
        w = 2 * intval(ones(1, dims(d)));
        w(1) = intval(1);
        shape = ones(1, dim);
        shape(d) = dims(d);
        nu_vect = nu_vect .* reshape(w, shape);
    end

    if dim >= 2
        % Last dimension (Chebyshev): weight = 2*nu^k (or 1 for k=0)
        cheb_idx = (0:dims(end)-1);
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
end

n = abs(a) .* nu_vect;

for i = 1:dim
    n = sum(n);
end

end
