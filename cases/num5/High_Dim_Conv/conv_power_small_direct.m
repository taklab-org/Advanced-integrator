function z = conv_power_small_direct(v, p)
% Direct convolution power for vectors and matrices.

if p == 2
    z = conv_small_direct(v, v);
elseif p == 3
    w = conv_small_direct(v, v);
    z = conv_small_direct(w, v);
else
    error('conv_power_small_direct: only p=2,3 supported');
end
end
