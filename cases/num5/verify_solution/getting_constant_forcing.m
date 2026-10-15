function [resolved,uncoupled] = getting_constant_forcing(f,m,h,params,L,nu,B,Vfin_w)
% Full-PDE Volterra bound for a time-independent, finitely supported forcing.
% f ALREADY includes its vector-field multiplier. Retain each tail mode's
% own integral; do not multiply the tail norm by its slowest decay rate.
f=intval(f);m=reshape(m,1,[]);n=prod(m);nf=size(f);
[kx,ky]=ndgrid(0:m(1)-1,0:m(2)-1);
w=2.^(double(kx>0)+double(ky>0)).*intval(nu).^(kx+ky);
fF=intval(zeros(m));overlap=min(m,nf);
fF(1:overlap(1),1:overlap(2))=f(1:overlap(1),1:overlap(2));
uncoupled=intval(zeros(n+1,1));
uncoupled(1:n)=intval(Vfin_w)*(w(:).*abs(fF(:)));
tail=f;tail(1:overlap(1),1:overlap(2))=0;
[tx,ty]=ndgrid(0:nf(1)-1,0:nf(2)-1);
outside=(tx>=m(1))|(ty>=m(2));
q=(intval(tx(outside))*intval(L(1))).^2+(intval(ty(outside))*intval(L(2))).^2;
d=intval(params(3))^2*q.^2-q+intval(params(2));
if ~all(inf(d)>0)
    error('getting_constant_forcing:Tail','Forcing tail is not strictly stable.');
end
hi=intval(h);
if isinf(sup(hi)),integral=1./d;else,integral=(1-exp(-d*hi))./d;end
tail(outside)=tail(outside).*integral;
uncoupled(n+1)=wnorm(tail,intval(nu));
resolved=intval(B)*uncoupled;
end
