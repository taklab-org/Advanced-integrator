function [psi] = psi_z_int(a,nu)
% [Opt5] Optimized: precompute iso_coeff, nu powers, move reshape outside loop.
N = size(a)-1;
n = N(end)+1;
numel_a = numel(mid(a));
psi = intval(zeros(numel_a ,1));
b = reshape(abs(a),[numel_a ,1]);

% [Opt5] Precompute all iso_coeff values
all_coeffs = zeros(numel_a, length(N));
for idx = 1:numel_a
    all_coeffs(idx, :) = iso_coeff(idx, N);
end
all_k = all_coeffs(:, 1:end-1);    % spatial indices
all_ell = all_coeffs(:, end);       % Chebyshev index

% [Opt5] Precompute nu powers for all possible ell2 values
max_ell2 = 2 * N(end);
nu_powers = intval(zeros(max_ell2 + 1, 1));
for p = 0:max_ell2
    nu_powers(p + 1) = nu^p;
end

for i = 2:numel_a
    k = all_k(i, :);
    ell = all_ell(i);
    for j = 1:numel_a
        j2 = all_k(j, :);
        ell2 = all_ell(j) + ell;
        if max(abs(k-j2)) <= N(1) || max(abs(k+j2)) <= N(1)
            if ell2 >= n && ell2 <= ell+n-1
                psi(j) = max(psi(j), b(j) / nu_powers(ell2 + 1));
            end
        end
    end
end
% [Opt5] reshape moved outside loop (was inside i-loop before)
psi = reshape(psi, size(a));
end
