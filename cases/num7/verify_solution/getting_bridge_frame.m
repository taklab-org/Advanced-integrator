function frame = getting_bridge_frame(abar,m,params,L,nu,kind)
% Fixed real basis in WEIGHTED finite coordinates; its 1-norm is <=1.
if nargin<6,kind='eigen';end
N=size(abar);v=abar(:);[~,A]=DF_OK_2D(v(2:end),params,N,L);
ids=reshape(1:prod(N),N);ids=ids(1:m(1),1:m(2));
[x,y]=ndgrid(0:m(1)-1,0:m(2)-1);
w=2.^(double(x>0)+double(y>0)).*intval(nu).^(x+y);D=diag(w(:));Di=diag(1./w(:));
if strcmp(kind,'physical'),V=eye(prod(m));
elseif strcmp(kind,'schur'),[V,~]=schur(mid(D*intval(A(ids(:),ids(:)))*Di),'real');
else,assert(strcmp(kind,'eigen'));[V,~]=eig(mid(D*intval(A(ids(:),ids(:)))*Di));end
if ~isreal(V),V=eye(prod(m));end
norms=sup(sum(abs(intval(V)),1));V=V./(norms*(1+1e-14));
V=intval(V);Vi=inv(V);
if ~all(isfinite(sup(abs(Vi))),'all'),V=intval(eye(prod(m)));Vi=V;end
assert(sup(norm(V,1))<=1);
frame=struct('V',V,'Vi',Vi,'D',D,'Di',Di,'m',m);
end
