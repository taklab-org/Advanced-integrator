function b = getting_block_norms(A,n)
% Weighted finite/tail operator block norms, with the tail stored as one norm.
A=abs(intval(A));b=intval(zeros(2));
b(1,1)=max(sum(A(1:n,1:n),1));b(1,2)=sum(A(1:n,n+1));
b(2,1)=max(A(n+1,1:n));b(2,2)=A(n+1,n+1);
b=intval(max(0,sup(b)));
end
