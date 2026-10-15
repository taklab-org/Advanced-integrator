function fx = eval_cheb(c,x,n)
% c: given (two-sided) chebyshev coefficients
% x: eval points in [a b]
arguments
    c; x; n = 200
end
a = x(1); b = x(end);
[M,m] = size(c);
k = 0:M-1;
xij = flipud(chebpts(n+1));
xc = (1.0 - xij)*a/2 + (1.0 + xij)*b/2; % Chebyshev points in [a,b]
fxc = cos(k.* acos(xij)) * c;
% interpolation via barycentric formula
valnum = length(x);
xi = 2*(x-a)/(b-a) - 1;
% xi = linspace(-1,1,valnum)';
% x = (1.0 - xi)*a/2 + (1.0 + xi)*b/2;
lam = [1/2; ones(n-1,1); 1/2] .* (-1).^((0:n)');
%
numer = zeros(valnum,m);
denom = zeros(valnum,1);
exact = zeros(valnum,1);
%
for j = 1:n+1
    xdiff = x - xc(j);
    temp = lam(j) ./ xdiff;
    numer = numer + temp * fxc(j,:);
    denom = denom + temp;
    exact(xdiff==0) = 1;
end
%
fx = numer ./ denom;
jj = find(exact);
fx(jj,:) = cos(k.* acos(xi(jj))) * c;
