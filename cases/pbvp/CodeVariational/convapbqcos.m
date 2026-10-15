function [vpshort,vpfull,vpextended] = convapbqcos(a,p,b,q)
N = size(a)-1;
dim=length(N);

p0=max(1,p); % takes care of p=0 case
q0=max(1,q); % takes care of p=0 case

for j=1:dim
  s{j}=1:N(j)+1;
  sflip{j}=abs(-N(j):N(j))+1;
  spfull{j}=(q0+p0)*N(j)+1:2*(q0+p0)*N(j)+1;
end

aa=a(sflip{1:dim});
bb=b(sflip{1:dim});

[~,vpextended]=convapbqfourier(aa,p,bb,q);
vpextended=real(vpextended);

vpfull = vpextended(spfull{1:dim});
vpshort = vpfull(s{1:dim});
end
