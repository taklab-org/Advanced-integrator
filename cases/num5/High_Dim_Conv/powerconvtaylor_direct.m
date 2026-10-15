function [vpshort,vpfull] = powerconvtaylor_direct(v,p)
% Direct convolution version of powerconvtaylor (no FFT).
% Direct mode is limited to vectors/matrices and p = 2 or 3.

if ndims(v) > 2
    error('powerconvtaylor_direct: only vectors and matrices are supported');
end

sz = size(v);
dim = length(sz);

vpfull = conv_power_small_direct(v, p);

for j = 1:dim
    s{j} = 1:sz(j);
end
vpshort = vpfull(s{1:dim});
end
