function [vpshort,vpfull,vpextended] = powerconvsin(v,p)
% computes the convolution power P for sine series V
% works in all dimensions
% Input:
%   v: sine coefficients, size (N1+1) x (N2+1) x ... (modes 0 to N per dimension)
%      Same convention as powerconvcos: first element is k=0 mode.
%      For pure sine, v(1,...) = 0 (no constant term).
%   p: convolution power (non-negative integer)
% Output (same format as powerconvcos):
%   vpshort: size of input (N+1 per dimension)
%   vpfull: positive half of full convolution, size (p*N+1) per dimension
%   vpextended: all Fourier modes, size (2*p*N+1) per dimension
%
% The antisymmetric extension is used: vv_k = sign(k)*v_{|k|}, vv_0=0.
%   For even p, the result is symmetric (cosine-type, k=0 nonzero).
%   For odd p, the result is antisymmetric (sine-type, k=0 is zero).
%   In both cases, the positive half [k=0, ..., k=pN] is returned as vpfull.
%
% Example (1D): v = [0, 1] represents sin(x)-type coefficient.
%   powerconvsin([0,1], 2) gives cosine coefficients of the squared result.

N = size(v)-1;
dim = length(N);

p0 = max(1,p); % takes care of p=0 case

% Build index arrays (same structure as powerconvcos)
for j=1:dim
  s{j} = 1:N(j)+1;                      % output short range (size N+1)
  sflip{j} = abs(-N(j):N(j))+1;         % symmetric folding index
  sfull{j} = p0*N(j)+1:2*p0*N(j)+1;     % positive half of extended
end

% Create symmetric extension (like powerconvcos)
vv = v(sflip{1:dim});

% Create sign tensor for antisymmetric extension
% sign_tensor has sign(k1)*sign(k2)*... for each dimension
% For N(j)=0 (trivial dimension), use 1 instead of sign(0)=0
sign_tensor = 1;
for j=1:dim
  krange = -N(j):N(j);
  if N(j) == 0
    sign_vals = 1;  % trivial dimension: no sign flip
  else
    sign_vals = sign(krange);
  end
  shape = ones(1, max(2,dim));  % at least 2 to avoid MATLAB squeezing
  shape(j) = 2*N(j)+1;
  sign_tensor = sign_tensor .* reshape(sign_vals, shape);
end

% Apply sign to get antisymmetric extension
if exist('intval','file') && isintval(v(1))
  % For intval: multiply sign tensor as intval
  sign_tensor = intval(sign_tensor);
end
vv = vv .* sign_tensor;

% Compute convolution power via Fourier
[~,vpextended] = powerconvfourier(vv,p);
if isreal(v)
  vpextended = real(vpextended);
end

% Extract positive half and short version
vpfull = vpextended(sfull{1:dim});
vpshort = vpfull(s{1:dim});

return
