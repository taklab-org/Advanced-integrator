function [vpshort,vpfull,vpextended,method] = powerconvsin_hybrid(v,p,opts)
% Hybrid powerconvsin choosing between FFT and direct convolution.

if nargin < 3
    opts = struct();
end

opts = powerconv_hybrid_defaults(opts);
method = powerconv_hybrid_method(v,p,opts);

switch method
    case 'direct'
        [vpshort,vpfull,vpextended] = powerconvsin_direct(v,p);
    case 'fft'
        [vpshort,vpfull,vpextended] = powerconvsin(v,p);
    otherwise
        error('powerconvsin_hybrid: unknown method %s', method);
end
end
