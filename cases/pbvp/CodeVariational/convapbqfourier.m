function [apbq,apbqfull] = convapbqfourier(a,p,b,q)
% computes the convolution power P of general Fourier series V
% works in all dimensions
% does not assume any symmetry
% the second, optional output is the "extended tensor"

sz = size(a);
dim=length(sz);

N = (sz-1)/2;

if p+q==0
    % special case of power=zero
    apbq=0*a;
    for j=1:dim
        s{j}=N(j)+1;
    end
    apbq(s{1:dim})=1;
    apbqfull=apbq;
    return
end
if p == 0
    a = 0*a;
    for j=1:dim
        s{j}=N(j)+1;
    end
    a(s{1:dim})=1;
end
if q == 0
    b = 0*b;
    for j=1:dim
        s{j}=N(j)+1;
    end
    b(s{1:dim})=1;
end

M = 2.^nextpow2(2*(p+q)*N+1);

a1 = altzeros(M(1:dim),a(1));
b1 = altzeros(M(1:dim),b(1));
M2=M/2;
M2(M==1)=0; % take care of dimensions with just one element (N=0)
s1 = M2 + 1 - N;
s2 = M2 + 1 + N;
s1f = M2 + 1 - (p+q)*N;
s2f = M2 + 1 + (p+q)*N;
for j=1:dim
  s{j}=s1(j):s2(j);
  sf{j}=s1f(j):s2f(j);
  fl{j}=[M2(j)+1:M(j),1:M2(j)];
end

a1(s{1:dim})=a;
a2 = a1(fl{1:dim}); %fftshift

b1(s{1:dim})=b;
b2 = b1(fl{1:dim}); %fftshift

clearvars -except a2 b2 p q sf s fl dim

alpha = altfftn(a2);
clearvars a2
beta = altfftn(b2);

w = (alpha.^p).*(beta.^q);


apbq2 = altifftn(w);

apbq1 = apbq2(fl{1:dim}); %fftshift

apbq = apbq1(s{1:dim});
apbqfull = apbq1(sf{1:dim});

return
