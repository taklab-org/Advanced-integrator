function [vpshort,vpfull,vpextended] = convapbqsin(a,p,b,q)
% convapbqsin  Convolution a^p * b^q for sine series.
%   Analogous to convapbqcos but with antisymmetric (sine) extension.
%   Both a and b are sine series coefficients [a_0; a_1; ...] (a_0=0).
%   For even p+q: result is cosine (symmetric) coefficients.
%   For odd p+q: result is sine (antisymmetric) coefficients.
%   Sign convention: same as powerconvsin (antisymmetric extension vv_k = sign(k)*v_{|k|}).

N = size(a)-1;
dim=length(N);

p0=max(1,p);
q0=max(1,q);

for j=1:dim
  s{j}=1:N(j)+1;
  sflip{j}=abs(-N(j):N(j))+1;
  spfull{j}=(q0+p0)*N(j)+1:2*(q0+p0)*N(j)+1;
end

% Symmetric folding (magnitude)
aa=a(sflip{1:dim});
bb=b(sflip{1:dim});

% Build sign tensor for antisymmetric extension
sign_tensor = 1;
for j=1:dim
  krange = -N(j):N(j);
  if N(j) == 0
    sign_vals = 1;
  else
    sign_vals = sign(krange);
  end
  shape = ones(1, max(2,dim));
  shape(j) = 2*N(j)+1;
  sign_tensor = sign_tensor .* reshape(sign_vals, shape);
end

if exist('intval','file') && (isintval(a(1)) || isintval(b(1)))
  sign_tensor = intval(sign_tensor);
end

% Apply antisymmetric extension
aa = aa .* sign_tensor;
bb = bb .* sign_tensor;

[~,vpextended]=convapbqfourier(aa,p,bb,q);
vpextended=real(vpextended);

vpfull = vpextended(spfull{1:dim});
vpshort = vpfull(s{1:dim});
end
