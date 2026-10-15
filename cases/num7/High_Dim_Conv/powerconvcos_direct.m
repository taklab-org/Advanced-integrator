function [vpshort,vpfull,vpextended] = powerconvcos_direct(v,p)
% Direct convolution version of powerconvcos (no FFT)
% Produces tighter interval bounds by avoiding frequency-domain amplification.
%
% v: 1D/2D one-sided cosine-Chebyshev coefficients
% p: power (2 or 3)
% vpfull: full output (p*n1-p+1 x p*n2-p+1)
% vpshort: truncated to input size

if ndims(v) > 2
    error('powerconvcos_direct: only vectors and matrices are supported');
end

N = size(v)-1;

% Create symmetric extension in both dimensions
sflip1 = abs(-N(1):N(1))+1;
sflip2 = abs(-N(2):N(2))+1;
vv = v(sflip1, sflip2); % size (2*N(1)+1, 2*N(2)+1)

z = conv_power_small_direct(vv, p);

if isreal(v)
    z = real(z);
end
vpextended = z;

% Extract positive-index portion (same as powerconvcos)
sfull1 = p*N(1)+1 : 2*p*N(1)+1;
sfull2 = p*N(2)+1 : 2*p*N(2)+1;
vpfull = vpextended(sfull1, sfull2);

s1 = 1:N(1)+1;
s2 = 1:N(2)+1;
vpshort = vpfull(s1, s2);
end
