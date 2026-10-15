function [vpshort,vpfull,vpextended] = powerconvcos(v,p)
% computes the convolution power P for cosine series V
% works in all dimensions
% vpfull is just the full cosine convolution
% vpextended is are all modes, not just the positive ones
% vpshort has size(v)
% if p=0, vpfull has size(v)

N = size(v)-1;
dim=length(N);

p0=max(1,p); % takes care of p=0 case

for j=1:dim
  s{j}=1:N(j)+1;
  sflip{j}=abs(-N(j):N(j))+1;
  sfull{j}=p0*N(j)+1:2*p0*N(j)+1;
end

vv=v(sflip{1:dim});

[~,vpextended]=powerconvfourier(vv,p);
vpextended=real(vpextended);

vpfull = vpextended(sfull{1:dim});
vpshort = vpfull(s{1:dim});

return
