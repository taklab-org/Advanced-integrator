function value = wopnorm(A,nu,q)% the weighted norm of matrix
arguments
  A, nu, q=0
end
[n,m] = size(A);
k = (0:n-1)';
ell = (0:m-1)';
alpk = 2*ones(n,1); alpk(1) = 1;
alpell = 2*ones(m,1); alpell(1) = 1;
w_k0 = alpell .* nu.^(ell);
w_kq = alpk .* nu.^(k) .* (1+k).^q;
%
value = max(sum(w_kq.*abs(A),1)./w_k0');
