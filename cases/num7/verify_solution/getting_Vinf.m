function [Vinf,Vinfq,mu_m] = getting_Vinf(h,N,m,params,L)
% Diagonal-tail Volterra bounds on the ENTIRE complement of F_m.
% A finite box is used only to sharpen Vinfq; its infinite complement is
% bounded analytically. h=Inf gives the infinite-time constants.
sigma=intval(params(2));eps2=intval(params(3))^2;
Li=intval(L);hi=intval(h);m=reshape(m,1,[]);N=reshape(N,1,[]);
assert(numel(m)==2 && numel(N)==2 && all(m>=1) && inf(hi)>=0);
r0=inf(min((intval(m(1))*Li(1))^2,(intval(m(2))*Li(2))^2));
r=intval(r0);d0=eps2*r^2-r+sigma;
% d(q)>0 and increasing; q/d(q) decreasing for every q>=r0.
if ~(inf(d0)>0 && inf(2*eps2*r-1)>0 && inf(eps2*r^2-sigma)>0)
    error('getting_Vinf:Tail','Infinite-tail monotonicity/stability could not be verified.');
end
mu_m=-d0; % upper bound for every tail eigenvalue, not a finite search.
if isinf(sup(hi))
    Vinf=1/d0;Vinfq=r/d0;return;
end
Vinf=(1-exp(-d0*hi))/d0;
if sup(hi)==0,Vinfq=intval(0);return;end
Nsearch=max(N,m+1);
for refinement=1:12
    [kx,ky]=ndgrid(0:Nsearch(1)-1,0:Nsearch(2)-1);
    outside=(kx>=m(1))|(ky>=m(2));
    q=(intval(kx(outside))*Li(1)).^2+(intval(ky(outside))*Li(2)).^2;
    d=eps2*q.^2-q+sigma;
    assert(all(inf(d)>0));
    finite_upper=sup(max(q.*(1-exp(-d*hi))./d));
    cut=inf(min((intval(Nsearch(1))*Li(1))^2,(intval(Nsearch(2))*Li(2))^2));
    qc=intval(cut);dc=eps2*qc^2-qc+sigma;
    % Beyond the searched rectangle, q>=cut and q/d(q)<=cut/d(cut).
    remainder_upper=sup(qc/dc);
    Vinfq=intval(max(finite_upper,remainder_upper));
    if remainder_upper<=finite_upper || prod(2*Nsearch)>2e6,return;end
    Nsearch=2*Nsearch;
end
% The last finite-box + analytic-complement bound is valid even if coarse.
end
