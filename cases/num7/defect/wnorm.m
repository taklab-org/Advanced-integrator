function norma = wnorm(a,nu,q) % return weighted ell-one norm
arguments
  a, nu, q=0;
end
% WeightedNORM
% a:  d-dimensional tensor (d: spatial dimension)
% nu: exponential weights (nu >= 1)
% q:  an algebraic weights for derivative of nonlinearity
%
if ismatrix(a)
    if isvector(a) % vector case
        n = length(a);
        k = (0:n-1)';
        alp = 2*ones(n,1); alp(1) = 1;
        w = alp .* nu.^(k) .* (1+k).^q;
    else  % matrix case
        [n1,n2] = size(a); % length of input
        [k1,k2] = ndgrid(0:n1-1,0:n2-1);
        alp = 4*ones(n1,n2); alp(1,:) = 2; alp(:,1) = 2; alp(1,1) = 1;
        v1 = nu.^(0:n1-1); v2 = nu.^(0:n2-1);
        w = alp .* (v1'*v2) .* (1+k1+k2).^q;
    end
else % 3D tensor case (TODO: more general code?)
    if exist('intval','file') && isintval(a(1))
        n = size(a); n1 = n(1); n2 = n(2); n3 = n(3);
    else
        [n1,n2,n3] = size(a); % length of input
    end
    [k1,k2,k3] = ndgrid(0:n1-1,0:n2-1,0:n3-1);
    alp = 8*ones(n1,n2,n3);
    alp(1,:,:) = 4; alp(:,1,:) = 4; alp(:,:,1) = 4;
    alp(1,1,:) = 2; alp(1,:,1) = 2; alp(:,1,1) = 2;
    alp(1,1,1) = 1;
    v1 = nu.^(0:n1-1); v2 = nu.^(0:n2-1); v3 = nu.^(0:n3-1);
    w = alp .* ((permute(v3,[3,1,2]).*ones(n1,n2,n3)) .* (v1'*v2)) .* (1+k1+k2+k3).^q;
end
%
norma = sum(abs(a).*w,"all");
