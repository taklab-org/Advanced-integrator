function [S,Cend,J] = getting_matrix_exp_bound(data,h,left,right)
% Signed midpoint Taylor bounds for left*V*exp(Lambda*t)*Vi*right.
% This is the SAME approximate exponential used by the Neumann correction.
n=size(data.V,1);if nargin<3,left=intval(eye(n));end
if nargin<4,right=intval(eye(n));end
P=intval(left)*data.V;R=data.Vi*intval(right);lam=data.lambda;
hi=intval(h);H=sup(hi);assert(isfinite(H) && H>=0);
Cend=P*diag(exp(lam*hi))*R;S=zeros(size(Cend));J=intval(zeros(size(Cend)));
rate=max(sup(abs(lam)));
if H==0,S=intval(sup(abs(intval(left)*intval(right))));return;end
if rate==0,nodes=[0,H];
else
    first=min(H,0.05/rate);count=max(0,ceil(log(H/first)/log(1.2)));
    nodes=unique([0,min(H,first*1.2.^(0:count)),H]);
end
for j=1:numel(nodes)-1
    ta=intval(nodes(j));tb=intval(nodes(j+1));tm=(ta+tb)/2;d=tb-ta;
    C=P*diag(exp(lam*tm))*R;Cp=P*diag(lam.*exp(lam*tm))*R;
    e=sup(exp(real(lam)*infsup(nodes(j),nodes(j+1))));
    C2=abs(P)*diag(abs(lam).^2.*intval(e))*abs(R);
    cellBound=abs(C)+(d/2).*abs(Cp)+(d^2/8).*C2; S=max(S,sup(cellBound)); if nargout>2,J=J+d*cellBound;end
end
S=intval(max(0,S));
end
