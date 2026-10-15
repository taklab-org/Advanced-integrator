function [r, Y_0, Z_0, Z_1] = eval_full_bounds(ctx, c, e_j)
Y_0 = eval_bound_y0(ctx, c, e_j);
Z_0 = ctx.Z_0;
Z_1 = ctx.Z_1;
den = 1 - Z_0 - Z_1;
if den <= 0
    error(strcat("Radii polynomial denominator 1 - Z_0 - Z_1 = ", num2str(den), " is not positive"))
end
r = Y_0 / den;
end
