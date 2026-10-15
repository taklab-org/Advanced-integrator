function [vpshort,vpfull,vpextended,method] = powerconvcos_hybrid(v,p,opts)
% Hybrid powerconvcos choosing between FFT and direct convolution.
% Direct mode is limited to vectors/matrices and p = 2 or 3.
%
% opts.mode: 'auto' (default), 'fft', or 'direct'
% opts.max_1d_p2: maximum vector length for direct p=2 (default 64)
% opts.max_1d_p3: maximum vector length for direct p=3 (default 32)
% opts.max_2d_double_p2: maximum matrix size for direct double p=2 (default 8)
% opts.max_2d_double_p3: maximum matrix size for direct double p=3 (default 8)
% opts.max_2d_interval_point_p2: maximum matrix size for point-interval direct p=2 (default 4)
% opts.max_2d_interval_point_p3: maximum matrix size for point-interval direct p=3 (default 4)
% opts.max_2d_interval_wide_p2: maximum matrix size for nonzero-radius interval direct p=2 (default 8)
% opts.max_2d_interval_wide_p3: maximum matrix size for nonzero-radius interval direct p=3 (default 6)

if nargin < 3
    opts = struct();
end

opts = powerconv_hybrid_defaults(opts);
method = powerconv_hybrid_method(v,p,opts);

switch method
    case 'direct'
        [vpshort,vpfull,vpextended] = powerconvcos_direct(v,p);
    case 'fft'
        [vpshort,vpfull,vpextended] = powerconvcos(v,p);
    otherwise
        error('powerconvcos_hybrid: unknown method %s', method);
end

end
