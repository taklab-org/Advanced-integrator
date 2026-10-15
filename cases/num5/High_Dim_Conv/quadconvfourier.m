function [ab,abfull] = quadconvfourier(a,b)
% computes the quadratic product of general Fourier series a,b
% works in all dimensions
% does not assume any symmetry
% the second, optional output is the "extended tensor"

if ~isequal(size(a), size(b))
    error('quadconvfourier: a and b must have the same size')
end

sz = size(a);
dim=length(sz);

if any(mod(sz,2)==0)
    error('quadconvfourier: input size must be odd in all dimensions (2*N+1)')
end

N = (sz-1)/2;

p = 2; % quadratic

M = 2.^nextpow2(2*p*N+1);

u = altzeros(M(1:dim),a(1));
v = altzeros(M(1:dim),b(1));

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

u(s{1:dim})=a; v(s{1:dim})=b;

u2 = u(fl{1:dim}); v2 = v(fl{1:dim}); %fftshift

uv = altfftn(u2).*altfftn(v2);

uv2 = altifftn(uv);

uv1 = uv2(fl{1:dim}); %fftshift

ab = uv1(s{1:dim});
abfull = uv1(sf{1:dim});

return
