function [DF] = DF_lin_fin_dim_KS(mu_k,Nop_coeff,a_bar,h,c,L,lambda)
% DF_lin_fin_dim_KS  Jacobian of variational equation for 1D KS (sine basis).
%   Branches 1-4 (Chebyshev differentiation, mu_k, boundary) are identical to SH.
%   Branch 5 (nonlinear DN contribution) uses sine convolution structure:
%     DN_{k1,k2} = k1*L * (a_bar_{k1+k2} + sgn(k2-k1)*a_bar_{|k1-k2|})
%   where a_bar is the sine Fourier coefficient (0-indexed: a_bar(1)=mode k=0=0).

if ~exist('lambda','var')
    lambda = [];
end

% --- Persistent cache ---
persistent cache
if local_cache_hit(cache, mu_k, Nop_coeff, a_bar, h, size(c), L, lambda)
    DF = cache.DF;
    return
end

N = size(a_bar)-1;
P = size(Nop_coeff,2)-1;
Mc = size(c);
Ma = size(a_bar);
dim = length(Mc);
dim_sp = dim - 1;
n_cheb = Mc(end);
num_sp = prod(Mc(1:end-1));
numel_c = prod(Mc);

% Preamble
mu_k_temp = mu_k;
mu_k_full = zeros(Ma);
for i = 1:length(N)
    otherdims_test{i} = 1:N(i)+1;
    otherdims_c{i} = 1:Mc(i);
end
mu_k_full(otherdims_c{:}) = mu_k_temp;

% Extended a_bar for spatial indexing (need up to 2*m modes)
N_ext_sp = 2*N(1);
a_bar_ext = zeros(N_ext_sp+1, N(end)+1);
a_bar_ext(1:N(1)+1, :) = a_bar;

% Strides for N-sized array
N_strides = ones(1, dim_sp);
for d = 2:dim_sp
    N_strides(d) = prod(N(1:d-1)+1);
end
cheb_stride_N = prod(N(1:dim_sp)+1);

% Spatial mode table (0-indexed wavenumbers)
sp_subs = cell(1, dim_sp);
sp_sz = Mc(1:end-1); if isscalar(sp_sz), sp_sz = [sp_sz, 1]; end
[sp_subs{:}] = ind2sub(sp_sz, (1:num_sp)');
k_sp = zeros(num_sp, dim_sp);
for d = 1:dim_sp
    k_sp(:,d) = sp_subs{d} - 1;  % 0-indexed: k=0,...,m-1
end

% Per spatial mode: K_q = k*L, mu_k base index
K_q_vals = zeros(num_sp, 1);
sp_mu_base = zeros(num_sp, 1);
for s = 1:num_sp
    K_q_vals(s) = k_sp(s,1) * L;  % k*L for 1D
    sp_mu_base(s) = 1 + sum(k_sp(s,:) .* N_strides);
end

% Lambda condition
skip_zero_sp = ~isempty(lambda);
is_zero_sp = all(k_sp == 0, 2);

% Chebyshev index matrices
[E1, E2] = ndgrid(0:n_cheb-1, 0:n_cheb-1);

% Nop_coeff factor for KS: i*Nop_coeff(j+1,i+1) combined with convapbqsin sign
% For Nop_coeff = [0,0,0; 0,0,-1/2]: factor = 2*(-1/2)*(-1) = 1
% where (-1) accounts for convapbqsin returning negative of true cosine coeffs.
% General: for each (i,j) entry, the Nop factor that multiplies the spatial DN entry.
nop_factor = 0;
for i = 1:(size(Nop_coeff,2)-1)
    for j = 0:(size(Nop_coeff,1)-1)
        if Nop_coeff(j+1,i+1) ~= 0
            % convapbqsin(a,i-1,c,1) returns (-1)^(i-1+1) * true result
            % For i-1+1 = i: even i gives +1 (true cosine), odd i gives -1
            % Actually: convapbqsin with total power i-1+1=i returns (-1)^i * true
            % Wait - the sign is (-1) for any sum of antisymmetric extensions squared
            % For p=i-1, q=1: total power = i. Sign = (-1)^(number of sine inputs) = (-1)^2 = 1 for even...
            % Actually: each antisymmetric ext introduces factor i in Fourier space.
            % Product of two antisymmetric exts: i^2 = -1. So output = -1 * true.
            % For convapbqsin(a, p, b, q) with total power p+q:
            %   each input contributes one i factor -> total (i)^(p+q) = (-1)^(p+q) for even p+q
            %   Wait no: it's i^2 = -1 (two sine inputs), regardless of powers.
            % convapbqsin(a, i-1, c, 1): both a and c get antisymmetric extension.
            % The output = (-1) * true cosine coefficients (for even total power).
            % So the effective factor is: i * Nop_coeff(j+1,i+1) * (-1)
            nop_factor = nop_factor + i * Nop_coeff(j+1,i+1) * (-1);
        end
    end
end

% Build DF block by block
DF = zeros(numel_c, numel_c);

for s1 = 1:num_sp
    idx_r = s1 + (0:n_cheb-1)*num_sp;
    active = ~skip_zero_sp || ~is_zero_sp(s1);

    % Branches 1-4: diagonal block (identical to SH)
    block_diag = zeros(n_cheb, n_cheb);
    block_diag(1,1) = 1;
    if n_cheb > 1
        block_diag(1, 2:end) = 2*(-1).^(1:n_cheb-1);
    end
    for ell1 = 1:n_cheb-1
        block_diag(ell1+1, ell1+1) = 2*ell1;
    end
    if active
        for ell1 = 1:n_cheb-1
            mu_val = mu_k_full(sp_mu_base(s1) + ell1*cheb_stride_N);
            block_diag(ell1+1, ell1) = block_diag(ell1+1, ell1) - h/2*mu_val;
            if ell1 <= n_cheb-2
                block_diag(ell1+1, ell1+2) = block_diag(ell1+1, ell1+2) + h/2*mu_val;
            end
        end
    end
    DF(idx_r, idx_r) = block_diag;

    % Branch 5: nonlinear (sine-specific DN structure)
    if active
        k1 = k_sp(s1, 1);  % 0-indexed spatial wavenumber
        K_q = K_q_vals(s1); % k1*L
        if K_q ~= 0  % skip k1=0 (K_q=0 → no nonlinear contribution)
            for s2 = 1:num_sp
                k2 = k_sp(s2, 1);  % 0-indexed
                if k2 == 0
                    continue;  % a_bar at k=0 is zero, no contribution
                end
                idx_c = s2 + (0:n_cheb-1)*num_sp;

                % Sine DN structure: DN_{k1,k2} = k1*L * (a_{k1+k2} + sgn(k2-k1)*a_{|k1-k2|})
                idx_sum = k1 + k2;     % Hankel index
                idx_diff = abs(k1 - k2); % Toeplitz index
                if k1 ~= k2
                    sign_diff = sign(k2 - k1);
                else
                    sign_diff = 0;  % ā_0 = 0 anyway
                end

                % Get Chebyshev series of DN coefficient
                % coeff_cheb = a_bar(idx_sum) + sign_diff * a_bar(idx_diff)
                % Multiply by nop_factor * K_q
                block_5 = zeros(n_cheb, n_cheb);
                for sigma3 = [+1, -1]
                    cheb_diff = E1 - sigma3*E2;

                    % Hankel part: a_bar at spatial index idx_sum
                    if idx_sum <= N_ext_sp
                        idx_p_h = abs(cheb_diff + 1);
                        idx_m_h = abs(cheb_diff - 1);
                        valid_h = (idx_p_h <= N(end)) & (idx_m_h <= N(end));
                        nop_p_h = zeros(n_cheb);
                        nop_m_h = zeros(n_cheb);
                        nop_p_h(valid_h) = a_bar_ext(idx_sum+1, idx_p_h(valid_h)+1);
                        nop_m_h(valid_h) = a_bar_ext(idx_sum+1, idx_m_h(valid_h)+1);
                        contrib_h = nop_factor * h/2 * K_q * (nop_p_h - nop_m_h);
                    else
                        contrib_h = zeros(n_cheb);
                    end

                    % Toeplitz part: sign_diff * a_bar at spatial index idx_diff
                    if sign_diff ~= 0 && idx_diff <= N_ext_sp
                        idx_p_t = abs(cheb_diff + 1);
                        idx_m_t = abs(cheb_diff - 1);
                        valid_t = (idx_p_t <= N(end)) & (idx_m_t <= N(end));
                        nop_p_t = zeros(n_cheb);
                        nop_m_t = zeros(n_cheb);
                        nop_p_t(valid_t) = a_bar_ext(idx_diff+1, idx_p_t(valid_t)+1);
                        nop_m_t(valid_t) = a_bar_ext(idx_diff+1, idx_m_t(valid_t)+1);
                        contrib_t = sign_diff * nop_factor * h/2 * K_q * (nop_p_t - nop_m_t);
                    else
                        contrib_t = zeros(n_cheb);
                    end

                    contrib = contrib_h + contrib_t;
                    contrib(1,:) = 0;
                    if sigma3 == -1
                        contrib(:,1) = 0;
                    end
                    block_5 = block_5 + contrib;
                end
                DF(idx_r, idx_c) = DF(idx_r, idx_c) + block_5;
            end
        end
    end
end

DF = reshape(DF, [Mc, Mc]);

% Store in cache
cache = struct();
cache.mu_k = mu_k_temp;
cache.Nop_coeff = Nop_coeff;
cache.a_bar = a_bar;
cache.h = h;
cache.c_size = size(c);
cache.L = L;
cache.lambda = lambda;
cache.DF = DF;
end

function hit = local_cache_hit(cache, mu_k, Nop_coeff, a_bar, h, c_size, L, lambda)
if isempty(cache)
    hit = false;
    return
end
hit = isequal(cache.c_size, c_size) ...
    && isequal(cache.h, h) ...
    && isequal(cache.L, L) ...
    && isequal(cache.lambda, lambda) ...
    && isequal(cache.mu_k, mu_k) ...
    && isequal(cache.Nop_coeff, Nop_coeff) ...
    && isequal(cache.a_bar, a_bar);
end
