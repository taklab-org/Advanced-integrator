function [vp,vpfull] = powerconvfourier(v,p)
% computes the convolution power P of general Fourier series V
% works in all dimensions
% does not assume any symmetry
% the second, optional output is the "extended tensor"

sz = size(v);
dim=length(sz);

if any(mod(sz,2)==0)
    error('powerconvfourier: input size must be odd in all dimensions (2*N+1)')
end

N = (sz-1)/2;

if p==0
    % special case of power=zero
    vp=0*v;
    for j=1:dim
        s{j}=N(j)+1;
    end
    vp(s{1:dim})=1;
    vpfull=vp;
    return
end

M = 2.^nextpow2(2*p*N+1);

v1 = altzeros(M(1:dim),v(1));
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

v1(s{1:dim})=v;
v2 = v1(fl{1:dim}); %fftshift

w = altfftn(v2);
wp = w.^p;
vp2 = altifftn(wp);

vp1 = vp2(fl{1:dim}); %fftshift

vp = vp1(s{1:dim});
vpfull = vp1(sf{1:dim});

return
