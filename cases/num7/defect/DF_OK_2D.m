function [DF,A] = DF_OK_2D(a,params,N,L)
% Jacobian of F_OK_2D, with the full square and all distinct source signs.
u=reshape([params(1);a],N);[~,b]=powerconvcos(u,2);
[x,y]=ndgrid(0:N(1)-1,0:N(2)-1);x=x(:);y=y(:);stride=size(b,1);
mm=abs(x-x.')+1+stride*abs(y-y.');
pm=x+x.'+1+stride*abs(y-y.');
mp=abs(x-x.')+1+stride*(y+y.');
pp=x+x.'+1+stride*(y+y.');
sx=double(x.'>0);sy=double(y.'>0);
DN=3*(b(mm)+b(pm).*sx+b(mp).*sy+b(pp).*(sx.*sy));
q=(x*L(1)).^2+(y*L(2)).^2;mu=q.*(1-params(3)^2*q)-params(2);
A=diag(mu)-q.*DN;DF=A(2:end,2:end);
end
