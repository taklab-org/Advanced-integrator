function Y_0 = eval_bound_y0(ctx, c, e_j)
% Compute only Y_0 using precomputed context + V2's vectorized norms.
% [Opt2] Uses precomputed nu_vect from ctx.
% [Opt3] Uses pre-allocated zero templates from ctx.
if ~isa(c, 'intval')
    c = intval(c);
end
if ~isa(e_j, 'intval')
    e_j = intval(e_j);
end

F_nm = reshape(F_lin_fin_dim_int(ctx.mu_k, ctx.Nop_coeff, ctx.a_bar, c, ctx.h, e_j, ctx.L, ctx.lambda), [ctx.numel_c, 1]);

% [Opt3] Reuse pre-allocated zero templates instead of intval(zeros(...))
c_ext = ctx.c_ext_zero;
c_ext(ctx.otherdims_c{:}) = c;
c_pad = ctx.c_pad_zero;
c_pad(ctx.otherdims_a{:}) = c_ext;
T_N_c_bar = T_op_int(ctx.h / 2 * N_op_lin_int(ctx.Nop_coeff, ctx.a_bar_pad, c_pad));

% [Opt2] Pass precomputed nu_vect as 3rd argument
Y_0_1 = norm_1_nu_int(reshape(ctx.A_nm * F_nm, ctx.Mc), ctx.nu, ctx.nu_vect_Mc);
lambda_k_c_k = ctx.mu_k .* c;
Y_0_2 = abs(ctx.nu ^ ctx.n / ctx.n * lambda_k_c_k(ctx.otherdims{:}, ctx.n));
for i = 1:length(ctx.N) - 1
    Y_0_2 = sum(Y_0_2);
end

Y_0_3 = intval(zeros(size(T_N_c_bar)));
for i = ctx.n:ctx.M
    Y_0_3(ctx.otherdims{:}, i + 1) = 1 / (2 * i) * T_N_c_bar(ctx.otherdims{:}, i + 1);
end
% [Opt2] Pass precomputed nu_vect for pad_size
Y_0_3 = norm_1_nu_int(Y_0_3, ctx.nu, ctx.nu_vect_pad);
Y_0 = sup(Y_0_1 + Y_0_2 + Y_0_3);
end
