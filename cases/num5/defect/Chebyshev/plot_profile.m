function plot_profile(c,I,n)
% c: given (two-sided) chebyshev coefficients (dim must be 2)
% I: domain of the function f
arguments
    c; I = [-1,1]; n = 200
end
a = I(1); b = I(2);
[M,m] = size(c);
if m ~= 2
    error('dim of data should be 2. If it is 3-dimensional, use plot3_profile instead.')
end
k = 0:M-1;
xij = flipud(chebpts(n+1));
xc = (1.0 - xij)*a/2 + (1.0 + xij)*b/2; % Chebyshev points in [a,b]
fxc = cos(k.* acos(xij)) * c;
% interpolation via barycentric formula
valnum = 5000;
xi = linspace(-1,1,valnum)';
x = (1.0 - xi)*a/2 + (1.0 + xi)*b/2;
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
%
plot(fx(:,1),fx(:,2),'linewidth',1.6)
xlabel('$x_1$','interpreter','latex')
ylabel('$x_2$','interpreter', 'latex')
% xlim(I)
