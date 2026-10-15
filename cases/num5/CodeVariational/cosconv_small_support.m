function result = cosconv_small_support(f, g, support_sz, out_sz)
%COSCONV_SMALL_SUPPORT  Direct 3D cosine convolution with small support.
%
%   result = cosconv_small_support(f, g, support_sz, out_sz)
%
%   Computes the cosine convolution (pointwise product in physical space)
%   of two 3D cosine coefficient arrays f and g, where g has nonzero
%   entries only at modes (0:m1-1, 0:m2-1, 0:ng-1).
%
%   This replaces the expensive 3D intval FFT in convapbqcos with a
%   direct summation that exploits the small spatial support of g.
%
%   Inputs:
%     f:          (F1, F2, Fn) intval array  (e.g., precomputed a_bar^2)
%     g:          (G1, G2, Gn) intval array  (e.g., c_pad, spatially sparse)
%     support_sz: [m1, m2, ng] actual nonzero support of g
%     out_sz:     [o1, o2, on] desired output size
%
%   The 3D cosine convolution formula:
%     h(k1,k2,ell) = sum_{j1,j2,p in signed support of g}
%       f(|k1-j1|, |k2-j2|, |ell-p|) * g(|j1|, |j2|, |p|)

m1 = support_sz(1);
m2 = support_sz(2);
ng = support_sz(3);

F1 = size(f, 1);
F2 = size(f, 2);
Fn = size(f, 3);

o1 = out_sz(1);
o2 = out_sz(2);
on = out_sz(3);

result = intval(zeros(o1, o2, on));

% Output Chebyshev indices (0-based)
ell = 0:on-1;

for j1 = 0:m1-1
    for j2 = 0:m2-1
        % Extract g's Chebyshev slice for this spatial mode
        g_cheb = reshape(g(j1+1, j2+1, 1:ng), 1, 1, []);  % (1,1,ng)

        % Skip if all zero
        g_sup = sup(abs(g_cheb));
        if all(g_sup(:) == 0), continue; end

        % Spatial sign combinations: all unique (+-j1, +-j2)
        signs = sym_dif([j1, j2]);

        for si = 1:size(signs, 1)
            s1 = signs(si, 1);
            s2 = signs(si, 2);

            % Source spatial indices: |k - s| + 1 (1-based) for output k = 0:o-1
            k1_src = abs((0:o1-1)' - s1) + 1;  % (o1, 1)
            k2_src = abs((0:o2-1)  - s2) + 1;   % (1, o2)

            % Validity masks
            v1 = k1_src <= F1;
            v2 = k2_src <= F2;
            k1s = min(k1_src, F1);
            k2s = min(k2_src, F2);

            % p = 0 term: h += g_0 * f(spatial, ell)
            idx0 = min(ell + 1, Fn);
            vld0 = (ell + 1) <= Fn;
            contrib = g_cheb(1) * f(k1s, k2s, idx0);
            if any(~vld0)
                contrib(:, :, ~vld0) = 0;
            end

            % p > 0 terms: h += g_p * [f(spatial, |ell-p|) + f(spatial, ell+p)]
            for p = 1:ng-1
                gp = g_cheb(p+1);
                if sup(abs(gp)) == 0, continue; end

                % |ell - p| indices (1-based into f's 3rd dim)
                idx_m = abs(ell - p) + 1;
                vld_m = idx_m <= Fn;
                idx_ms = min(idx_m, Fn);

                % ell + p indices (1-based)
                idx_p = ell + p + 1;
                vld_p = idx_p <= Fn;
                idx_ps = min(idx_p, Fn);

                f_m = f(k1s, k2s, idx_ms);
                if any(~vld_m)
                    f_m(:, :, ~vld_m) = 0;
                end

                f_p = f(k1s, k2s, idx_ps);
                if any(~vld_p)
                    f_p(:, :, ~vld_p) = 0;
                end

                contrib = contrib + gp * (f_m + f_p);
            end

            % Zero out invalid spatial modes
            if any(~v1), contrib(~v1, :, :) = 0; end
            if any(~v2), contrib(:, ~v2, :) = 0; end

            result = result + contrib;
        end
    end
end
end
