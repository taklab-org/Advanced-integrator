function [vpshort,vpfull,vpextended] = powerconvsin_direct(v,p)
% Direct convolution version of powerconvsin (no FFT).
% Direct mode is limited to vectors/matrices and p = 2 or 3.

if ndims(v) > 2
    error('powerconvsin_direct: only vectors and matrices are supported');
end

N = size(v)-1;
dim = length(N);

for j = 1:dim
    s{j} = 1:N(j)+1;
    sflip{j} = abs(-N(j):N(j)) + 1;
    sfull{j} = p*N(j)+1 : 2*p*N(j)+1;
end

vv = v(sflip{1:dim});

sign_tensor = 1;
for j = 1:dim
    krange = -N(j):N(j);
    if N(j) == 0
        sign_vals = 1;
    else
        sign_vals = sign(krange);
    end
    shape = ones(1, max(2, dim));
    shape(j) = 2*N(j)+1;
    sign_tensor = sign_tensor .* reshape(sign_vals, shape);
end

if exist('intval','file') && isintval(v(1))
    sign_tensor = intval(sign_tensor);
end

vv = vv .* sign_tensor;
vpextended = conv_power_small_direct(vv, p);

if isreal(v)
    vpextended = real(vpextended);
end

vpfull = vpextended(sfull{1:dim});
vpshort = vpfull(s{1:dim});
end
