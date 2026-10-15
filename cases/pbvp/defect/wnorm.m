function norma = wnorm(a,nu,q) % return weighted ell-one norm
arguments
  a, nu, q=0;
end
% WeightedNORM
% a:  d-dimensional tensor (d: spatial dimension)
% nu: exponential weights (nu >= 1)
% q:  an algebraic weights for derivative of nonlinearity
%
inu = intval(nu);
if ismatrix(a)
    if isvector(a) % vector case
        n = length(a);
        k = (0:n-1)';
        alp = 2*ones(n,1); alp(1) = 1;
        w = intval(alp) .* inu.^(k) .* intval((1+k).^q);
    else  % matrix case
        [n1,n2] = size(a); % length of input
        [k1,k2] = ndgrid(0:n1-1,0:n2-1);
        alp = 4*ones(n1,n2); alp(1,:) = 2; alp(:,1) = 2; alp(1,1) = 1;
        v1 = inu.^(0:n1-1); v2 = inu.^(0:n2-1);
        w = intval(alp) .* (v1'*v2) .* intval((1+k1+k2).^q);
    end
else % 3D tensor case (TODO: more general code?)
    n = size(a); % length of input
    n1 = n(1); n2 = n(2); n3 = n(3);
    [k1,k2,k3] = ndgrid(0:n1-1,0:n2-1,0:n3-1);
    alp = 8*ones(n1,n2,n3);
    alp(1,:,:) = 4; alp(:,1,:) = 4; alp(:,:,1) = 4;
    alp(1,1,:) = 2; alp(1,:,1) = 2; alp(:,1,1) = 2;
    alp(1,1,1) = 1;
    v1 = inu.^(0:n1-1); v2 = inu.^(0:n2-1); v3 = inu.^(0:n3-1);
    w = intval(alp) .* ((permute(v3,[3,1,2]).*ones(n1,n2,n3)) .* (v1'*v2)) .* intval((1+k1+k2+k3).^q);
end
%
norma = sum(abs(a).*w,"all");
