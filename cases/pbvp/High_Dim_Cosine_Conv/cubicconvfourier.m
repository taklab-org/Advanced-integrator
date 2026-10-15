function [abc,abcfull] = cubicconvfourier(a,b,c)
% computes the cubic product of general Fourier series a,b,c
% works in all dimensions
% does not assume any symmetry
% the second, optional output is the "extended tensor"

sz = size(a);
dim=length(sz);

N = (sz-1)/2;

p = 3; % cubic

M = 2.^nextpow2(2*p*N+1);

u = altzeros(M(1:dim),a(1));
v = altzeros(M(1:dim),b(1));
w = altzeros(M(1:dim),c(1));

M2=M/2;
M2(M==1)=0; % take care of dimensions with just one element (N=0)
s1 = M2 + 1 - N;
s2 = M2 + 1 + N;
s1f = M2 + 1 - p*N;
s2f = M2 + 1 + p*N;
for j=1:dim
  s{j}=s1(j):s2(j);
  sf{j}=s1f(j):s2f(j);
  fl{j}=[M2(j)+1:M(j),1:M2(j)];
end

u(s{1:dim})=a; v(s{1:dim})=b; w(s{1:dim})=c;

u2 = u(fl{1:dim}); v2 = v(fl{1:dim});  w2 = w(fl{1:dim}); %fftshift

uvw = altfftn(u2).*altfftn(v2).*altfftn(w2);

uvw2 = altifftn(uvw);

uvw1 = uvw2(fl{1:dim}); %fftshift

abc = uvw1(s{1:dim});
abcfull = uvw1(sf{1:dim});

return
