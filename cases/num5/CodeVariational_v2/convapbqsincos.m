function [vpshort,vpfull,vpextended] = convapbqsincos(a,p,b,q)
% convapbqsincos  Convolution a^p * b^q for sine-Chebyshev series.
%   Spatial dimensions (all except last) use antisymmetric (sine) extension.
%   Last dimension (Chebyshev/temporal) uses symmetric (cosine) extension.
%   For 1D input, falls back to convapbqsin behavior (all sine).

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

% Build sign tensor: antisymmetric for spatial dims (1:dim-1), symmetric for last dim
if dim >= 2
  sign_tensor = 1;
  for j=1:dim-1
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

  aa = aa .* sign_tensor;
  bb = bb .* sign_tensor;
else
  % 1D: pure sine (same as convapbqsin)
  krange = -N(1):N(1);
  if N(1) == 0
    sign_vals = 1;
  else
    sign_vals = sign(krange);
  end

  if exist('intval','file') && (isintval(a(1)) || isintval(b(1)))
    sign_vals = intval(sign_vals);
  end

  aa = aa .* sign_vals(:);
  bb = bb .* sign_vals(:);
end

[~,vpextended]=convapbqfourier(aa,p,bb,q);
vpextended=real(vpextended);

vpfull = vpextended(spfull{1:dim});
vpshort = vpfull(s{1:dim});
end
