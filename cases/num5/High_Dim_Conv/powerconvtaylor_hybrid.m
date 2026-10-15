function [vpshort,vpfull,method] = powerconvtaylor_hybrid(v,p,opts)
% Hybrid powerconvtaylor choosing between FFT and direct convolution.

if nargin < 3
    opts = struct();
end

opts = powerconv_hybrid_defaults(opts);
method = powerconv_hybrid_method(v,p,opts);

switch method
    case 'direct'
        [vpshort,vpfull] = powerconvtaylor_direct(v,p);
    case 'fft'
        [vpshort,vpfull] = powerconvtaylor(v,p);
    otherwise
        error('powerconvtaylor_hybrid: unknown method %s', method);
end
end
