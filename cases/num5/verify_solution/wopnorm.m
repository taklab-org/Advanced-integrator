function value = wopnorm(A,nu,q)% the weighted norm of matrix
arguments
  A, nu, q=0
end
n = size(A,1);
k = (0:n-1)';
alp = 2*ones(n,1); alp(1) = 1;
w_k0 = alp .* nu.^(k);
w_kq = alp .* nu.^(k) .* (1+k).^q;
%
value = max(sum(w_kq.*abs(A),1)./w_k0');
