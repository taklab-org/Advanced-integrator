function c = conv_small_direct(a, b)
% Direct full convolution for vectors and matrices.

if ndims(a) > 2 || ndims(b) > 2
    error('conv_small_direct: only vectors and matrices are supported');
end

[m1, m2] = size(a);
[n1, n2] = size(b);
p1 = m1 + n1 - 1;
p2 = m2 + n2 - 1;

if exist('intval','file') && (isintval(a(1)) || isintval(b(1)))
    c = intval(zeros(p1, p2));
else
    c = zeros(p1, p2);
end

for i = 1:m1
    for j = 1:m2
        c(i:i+n1-1, j:j+n2-1) = c(i:i+n1-1, j:j+n2-1) + a(i,j) * b;
    end
end
end
