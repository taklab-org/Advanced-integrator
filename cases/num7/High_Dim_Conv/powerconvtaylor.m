function [vpshort,vpfull] = powerconvtaylor(v,p)
% computes the convolution power P for Taylor/power series V
% works in all dimensions
% Input:
%   v: Taylor coefficients, size (N1+1) x (N2+1) x ... (modes 0 to N per dimension)
%   p: convolution power (non-negative integer)
% Output:
%   vpshort: same size as v (N+1 per dimension)
%   vpfull: full convolution result (p*N+1 per dimension)
%
% Example: v = [1, 1] represents 1 + x
%   powerconvtaylor(v, 2) gives [1, 2, 1] representing (1+x)^2 = 1 + 2x + x^2
%
% Note: Taylor series are one-sided (k >= 0), so no fftshift is needed.
%   The k=0 coefficient is at index 1 in MATLAB convention, which aligns
%   with FFT's zero-frequency position.

sz = size(v);
dim = length(sz);

N = sz - 1;

if p==0
    % special case: power=zero returns 1 at k=0
    vp = 0*v;
    for j=1:dim
        s0{j} = 1;
    end
    vp(s0{1:dim}) = 1;
    vpshort = vp;
    vpfull = vp;
    return
end

% Pad size: need at least p*N+1 per dimension, round up to power of 2
M = 2.^nextpow2(p*N+1);
% Handle dimensions where N=0 (singleton): M=1 is fine (2^0=1)

v1 = altzeros(M(1:dim), v(1));

for j=1:dim
    s{j} = 1:sz(j);       % original range (1 to N+1)
    sf{j} = 1:p*N(j)+1;   % full output range (1 to p*N+1)
end

% Embed Taylor coefficients at the beginning (no shift needed)
v1(s{1:dim}) = v;

% FFT-based convolution: FFT -> power -> IFFT
w = altfftn(v1);
wp = w.^p;
vp1 = altifftn(wp);

vpfull = vp1(sf{1:dim});
vpshort = vp1(s{1:dim});

if isreal(v)
    vpfull = real(vpfull);
    vpshort = real(vpshort);
end

return
