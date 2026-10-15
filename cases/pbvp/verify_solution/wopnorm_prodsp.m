function value_mat = wopnorm_prodsp(A,nu,Udims,dif)% the weighted norm of matrix
arguments
  A, nu, Udims,dif=0 % default is ||G||
end

if dif
  % opnorm of B
  A11 = A(1:Udims,1); A12 = A(1:Udims,2:end);
  A21 = A(Udims+1:end,1);   A22 = A(Udims+1:end,2:end);
  value_mat = [norm(A11,1), wopnorm(A12,nu); ...
    wnorm(A21,nu), max(wopnorm(A22,nu),1)];
else
  % opnorm of G
  A11 = A(1:Udims,1:Udims); A12 = A(1:Udims,Udims+1:end);
  A21 = A(Udims+1:end,1:Udims); A22 = A(Udims+1:end,Udims+1:end);
  value_mat = [norm(A11,1), wopnorm(A12,nu); ...
    max([wnorm(A21(:,1),nu), wnorm(A21(:,2),nu), wnorm(A21(:,3),nu), wnorm(A21(:,4),nu)]), wopnorm(A22,nu)];
end
