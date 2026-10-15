function [abcshort,abcfull,abcextended] = cubicconvcos(a,b,c)
% computes the convolution power P for cosine series a,b,c
% works in all dimensions
% abcfull is just the full cosine convolution
% abcextended is are all modes, not just the positive ones
% abcshort has size(a)

N = size(a)-1;
dim=length(N);

p = 3; % cubic
p0=max(1,p); % takes care of p=0 case

for j=1:dim
  s{j}=1:N(j)+1;
  sflip{j}=abs(-N(j):N(j))+1;
  sfull{j}=p0*N(j)+1:2*p0*N(j)+1;
end

aa=a(sflip{1:dim});
bb=b(sflip{1:dim});
cc=c(sflip{1:dim});

[~,abcextended]=cubicconvfourier(aa,bb,cc);
abcextended=real(abcextended);

abcfull = abcextended(sfull{1:dim});
abcshort = abcfull(s{1:dim});

return
