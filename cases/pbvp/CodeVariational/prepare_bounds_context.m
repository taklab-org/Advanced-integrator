function ctx = prepare_bounds_context(mu_k, Nop_coeff, a_bar, h, nu, DF_nm, L, lambda)
% Precompute Z_0, Z_1, A_nm once (from New1), using V2's vectorized norms.
if ~exist('lambda', 'var')
    lambda = [];
end

N = size(mu_k) - 1;
n = N(end) + 1;
M = (size(Nop_coeff, 2) - 1) * (n - 1) + 1;
Mc = size(mu_k);
Ma = size(a_bar);
pad_size = [Ma(1:end-1), M + 1];

otherdims = repmat({':'}, 1, length(N) - 1);
otherdims_a = cell(1, length(Ma));
otherdims_c = cell(1, length(Mc));
for i = 1:length(Ma)
    otherdims_a{i} = 1:Ma(i);
    otherdims_c{i} = 1:Mc(i);
end

a_bar_pad = intval(zeros(pad_size));
a_bar_pad(otherdims_a{:}) = a_bar;

A_nm = inv(mid(DF_nm));
B = reshape(eye(size(A_nm)) - A_nm * DF_nm, [N + 1, N + 1]);
Z_0 = sup(norm_B_nu_int(B, nu, N));
Z_1 = local_compute_z1(mu_k, Nop_coeff, a_bar, h, nu, A_nm, L, N, n, otherdims, otherdims_c);
if Z_1 > 1
    error(strcat("Bound Z_1 = ", num2str(Z_1), " is > 1"))
end

ctx = struct();
ctx.mu_k = mu_k;
ctx.Nop_coeff = Nop_coeff;
ctx.a_bar = a_bar;
ctx.h = h;
ctx.nu = nu;
ctx.DF_nm = DF_nm;
ctx.L = L;
ctx.lambda = lambda;
ctx.N = N;
ctx.n = n;
ctx.M = M;
ctx.Mc = Mc;
ctx.Ma = Ma;
ctx.numel_c = numel(mid(mu_k));
ctx.otherdims = otherdims;
ctx.otherdims_a = otherdims_a;
ctx.otherdims_c = otherdims_c;
ctx.pad_size = pad_size;
ctx.a_bar_pad = a_bar_pad;
ctx.A_nm = A_nm;
ctx.Z_0 = Z_0;
ctx.Z_1 = Z_1;

% [Opt2] Precompute nu_vect for two array sizes used in eval_bound_y0
ctx.nu_vect_Mc = local_build_nu_vect(nu, Mc);
ctx.nu_vect_pad = local_build_nu_vect(nu, pad_size);

% [Opt3] Pre-allocate zero templates for eval_bound_y0
ctx.c_ext_zero = intval(zeros(Ma));
ctx.c_pad_zero = intval(zeros(pad_size));
end

function nv = local_build_nu_vect(nu, dims)
dim = length(dims);
w1 = 2 * intval(ones(dims(1), 1));
w1(1) = intval(1);
nv = w1;
for d = 2:dim-1
    w = 2 * intval(ones(1, dims(d)));
    w(1) = intval(1);
    shape = ones(1, dim);
    shape(d) = dims(d);
    nv = nv .* reshape(w, shape);
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
    nv = nv .* reshape(w_cheb, shape);
end
end

function Z_1 = local_compute_z1(mu_k, Nop_coeff, a_bar, h, nu, A_nm, L, N, n, otherdims, otherdims_c)
a_norm = norm_1_nu_int(a_bar, nu);
Z_1_2 = 1 / (2 * n) * (nu + 1 / nu) * (h / 2 * abs(mu_k(end)));
for i = 1:size(Nop_coeff, 2) - 1
    for j = 0:size(Nop_coeff, 1) - 1
        if Nop_coeff(j + 1, i + 1) ~= 0
            if j == 0
                Z_1_2 = Z_1_2 + 1 / (2 * n) * (nu + 1 / nu) * (h * abs(Nop_coeff(i + 1)) * a_norm ^ (i - 1));
            else
                NN = size(a_bar) - 1;
                dim = length(NN) - 1;
                s = arrayfun(@(k) 0:k, NN(1:end-1), 'UniformOutput', false);
                k = cell(1, dim);
                [k{:}] = ndgrid(s{:});
                K_sq = zeros(size(k{1}));
                for ell = 1:dim
                    K_sq = K_sq + (k{ell} * L(ell)) .^ j;
                end
                norm_ka = norm_1_nu_int(K_sq .* a_bar, nu);
                Z_1_2 = Z_1_2 + 1 / (2 * n) * (nu + 1 / nu) * (h * abs(Nop_coeff(i + 1)) * norm_ka ^ (i - 1));
            end
        end
    end
end

a_bar_z1 = a_bar;
a_bar_z1(otherdims{:}, N(end) + 2) = 0;
one_vect = zeros(size(a_bar_z1));
one_vect(1) = 1;
nop = N_op_lin_int(Nop_coeff, a_bar_z1, one_vect);
nop = nop(otherdims_c{:});

N_nop = size(nop) - 1;
dim = length(N_nop) - 1;
s = arrayfun(@(k) 0:k, N_nop(1:end-1), 'UniformOutput', false);
k = cell(1, dim);
[k{:}] = ndgrid(s{:});
K_sq = zeros(size(k{1}));
for ell = 1:dim
    K_sq = K_sq + (k{ell} * L(ell)) .^ (size(Nop_coeff, 1) - 1);
end
N_rep = ones(size(N_nop));
N_rep(end) = N_nop(end) + 1;
K_sq = repmat(K_sq, N_rep);
nop = K_sq .* nop;

psi = psi_z_int(nop, nu);
z = h * T_op_abs_int(psi);
z(otherdims{:}, 1) = 1 / (nu ^ n);
z = reshape(z(otherdims_c{:}), [prod(N + 1), 1]);
Z_1_1 = norm_1_nu_int(reshape(abs(A_nm) * z, N + 1), nu);
Z_1 = sup(Z_1_1 + Z_1_2);
end
