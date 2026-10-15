function [abshort,abfull,abextended] = quadconvcos(a,b)
% computes the quad convolution for cosine series a,b
% works in all dimensions
% abfull is just the full cosine convolution
% abextended is are all modes, not just the positive ones
% abshort has size(a)

N = size(a)-1;
dim=length(N);

p = 2; % quad
p0=max(1,p); % takes care of p=0 case

for j=1:dim
  s{j}=1:N(j)+1;
  sflip{j}=abs(-N(j):N(j))+1;
  sfull{j}=p0*N(j)+1:2*p0*N(j)+1;
end

aa=a(sflip{1:dim});
bb=b(sflip{1:dim});
% cc=c(sflip{1:dim});

[~,abextended]=quadconvfourier(aa,bb);
abextended=real(abextended);

abfull = abextended(sfull{1:dim});
abshort = abfull(s{1:dim});

return
