function fvals = chebvals(c,x,I)
% Evaluate two-sided Chebyshev coefficients by Clenshaw recurrence.
% This avoids acos of outward-rounded intervals at domain endpoints.
arguments
    c, x
    I = [-1,1];
end
if isempty(x),fvals=[];return;end
isiv=isa(c,'intval') || isa(x,'intval') || isa(I,'intval');
if isiv,c=intval(c);x=intval(x);I=intval(I);end
xi=2*(x(:)-I(1))/(I(2)-I(1))-1;
shape=[size(xi,1),size(c,2)];
if isiv,b1=intval(zeros(shape));else,b1=zeros(shape);end
b2=b1;
for j=size(c,1):-1:2
    b0=2*xi.*b1-b2+c(j,:);b2=b1;b1=b0;
end
if isempty(c),fvals=b1;else,fvals=xi.*b1-b2+c(1,:);end
end
