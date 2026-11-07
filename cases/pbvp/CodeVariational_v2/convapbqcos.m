function [vpshort,vpfull,vpextended] = convapbqcos(a,p,b,q)
% Same finite cosine convolution, with FFT padding from exact support only.
% No coefficient is discarded unless its absolute interval upper bound is 0.
% Supports are symmetric after cosine extension. Their Minkowski sum has
% half-width p*Ra+q*Rb, so M>=2*(p*Ra+q*Rb)+1 prevents circular aliasing.
% All nonzero interval calculations use the original verified FFT kernels.
N=size(a)-1; dim=numel(N);
assert(isequal(size(a),size(b)) && p>=0 && q>=0 && p==fix(p) && q==fix(q));
if p+q==0
    sflip=arrayfun(@(n) abs(-n:n)+1,N,'UniformOutput',false);
    spfull=arrayfun(@(n) 2*n+1:4*n+1,N,'UniformOutput',false);
    [~,vpextended]=convapbqfourier(a(sflip{:}),p,b(sflip{:}),q);
    vpextended=real(vpextended);vpfull=vpextended(spfull{:});
    s=arrayfun(@(n) 1:n+1,N,'UniformOutput',false);vpshort=vpfull(s{:});return
end
[Ra,zeroA]=exact_support(a,p,dim);[Rb,zeroB]=exact_support(b,q,dim);
R=p*Ra+q*Rb; M=2.^nextpow2(2*R+1);
type_sample=a(1)+b(1);
if zeroA || zeroB
    small=altzeros(R+1,type_sample);
else
    if p>0,alpha=cosine_fft(a,Ra,M,dim);else,alpha=1;end
    if q>0,beta=cosine_fft(b,Rb,M,dim);else,beta=1;end
    values=altifftn((alpha.^p).*(beta.^q));
    positive=arrayfun(@(r) 1:r+1,R,'UniformOutput',false);
    small=real(values(positive{:}));
end
vpshort=altzeros(N+1,type_sample);
overlap=min(N,R);s=arrayfun(@(r) 1:r+1,overlap,'UniformOutput',false);
vpshort(s{:})=small(s{:});
if nargout>=2
    F=(max(1,p)+max(1,q))*N;vpfull=altzeros(F+1,type_sample);
    overlap=min(F,R);s=arrayfun(@(r) 1:r+1,overlap,'UniformOutput',false);
    vpfull(s{:})=small(s{:});
end
if nargout>=3
    F=(p+q)*N;vpextended=altzeros(2*F+1,type_sample);
    target=arrayfun(@(f,r) f+1-r:f+1+r,F,R,'UniformOutput',false);
    mirrored=arrayfun(@(r) abs(-r:r)+1,R,'UniformOutput',false);
    vpextended(target{:})=small(mirrored{:});
end
end

function [R,iszero]=exact_support(a,p,dim)
R=zeros(1,dim);iszero=false;
if p==0,return;end
if isa(a,'intval'),bounds=sup(abs(a));else,bounds=abs(a);end
occupied=find(bounds>0);
if isempty(occupied),iszero=true;return;end
subs=cell(1,dim);[subs{:}]=ind2sub(size(a),occupied);
for d=1:dim,R(d)=max(subs{d})-1;end
end

function alpha=cosine_fft(a,R,M,dim)
mirrored=arrayfun(@(r) abs(-r:r)+1,R,'UniformOutput',false);
target=arrayfun(@(r,m) mod(-r:r,m)+1,R,M,'UniformOutput',false);
work=altzeros(M,a(1));work(target{:})=a(mirrored{:});
alpha=altfftn(work);
end
