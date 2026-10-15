function [Ta,Ma] = multiplication(a,M,Nout)
% [TENSOR,MATRIX] = multiplication(V,M,N)
% computes the tensor TENSOR and MATRIX of size N
% representing (left) multiplication by V,
% which is of Fourier-Chebyshev-Taylor type as encoded by M.
%
% Works for both floats and intervals
%
% If ndims(V)==sum(M)+2 then the first two dimensions are assumed
% to store the matrix dimensions (LxK) of a Jacobian of a vector field.
%
% N should be a matrix with two rows so that
% TENSOR is a tensor of size L x K x N(1,:) x N(2,:)
% and MATRIX is a matrix of size L*N(1,:) x K*N(2,:)
% where N(i,:) should be interpreted appropriately in both lines.
% The first row of N determines the "number of rows "in TENSOR/MATRIX
% The second row of N determines the "number of colums" in TENSOR/MATRIX
% If N consists of a single row it is assumed that N(2,:)=N(1,:)
% and TENSOR/MATRIX is square for K=L.

% setting some sizes
% square matrix case
if size(Nout,1)==1
  Nout(2,:)=Nout(1,:);
end
M=reshape(M,1,[]);
if length(M)<3
  M = [M zeros(1,3-length(M))];
end
sizeM = zeros(size(Nout,1),M(1)+M(2)+M(3));
for m=1:M(1)
  sizeM(:,m) = 2*Nout(:,m)+1;
end
for m=M(1)+1:M(1)+M(2)+M(3)
  sizeM(:,m) = Nout(:,m)+1;
end

% multiplication by scalar Fourier-Chebyshev-Taylor
% as special case of matrix-valued multiplation
if length(size(a))==sum(M)
  a=reshape(a,[1,1,size(a)]);
elseif length(size(a))==sum(M)+1
  if sum(M)==1 && size(a,2)==1
    % just a column vector of coefficients
    a=reshape(a,[1,1,size(a)]);
  else
    a=reshape(a,[1,size(a)]);
  end
end

sizea=size(a);
a=reshape(a,[sizea(1:2),prod(sizea(3:end))]);
prodsizea=[cumprod(sizea(4:end),"reverse"),1];

% preallocate
Ta = altzeros([sizea(1),sizea(2),prod(sizeM(1,:)),prod(sizeM(2,:))],a(1));
intervals = (exist('intval','file') && isintval(a(1)));

for d=1:sizea(1)
  for dd=1:sizea(2)
    % each derivative separately
    Mad=a(d,dd,:);
    sizeMad=[1;1];
    for m=1:M(1)
      % we do this recursively
      % Fourier
      Na=(sizea(m+2)-1)/2;
      K=(sizeM(:,m)-1)/2;
      NK1=min(Na,2*K(1));
      NK2=min(Na,2*K(2));
      % toeplitz matrix of indices
      T_K = toeplitz([(Na+1:Na+1+NK1) zeros(1,2*K(1)-NK1)]',[(Na+1:-1:Na+1-NK2) zeros(1,2*K(2)-NK2)]);
      % reshape with next dimension in third dimension of the tensor
      Mad = reshape(Mad,[sizeMad(1:2)',sizea(2+m),prodsizea(m)]);
      % concatenate zero matrix
      Madzero = altzeros([size(Mad,1:2),1,size(Mad,4)],a(1));
      % cat(3,A,B) did not always work nicely for intervals for with singleton dimensions
      Mad = permute([permute(Madzero,[3 1 2 4]);permute(Mad,[3 1 2 4])],[2 3 1 4]);
      % Mad = cat(3,Madzero,Mad);
      if ~intervals
        % turn into cell array of matrices
        Mad = mat2cell(Mad, sizeMad(1), sizeMad(2), ones(2*Na+2,1), size(Mad,4));
        % take toeplitz of this array and convert back to cells
        Mad = cell2mat(Mad(T_K+1));
      else
        % for intervals split into midpoint and radius
        Mad_mid = mat2cell(mid(Mad), sizeMad(1), sizeMad(2), ones(2*Na+2,1), size(Mad,4));
        Mad_mid = cell2mat(Mad_mid(T_K+1));
        Mad_rad = mat2cell(rad(Mad), sizeMad(1), sizeMad(2), ones(2*Na+2,1), size(Mad,4));
        Mad_rad = cell2mat(Mad_rad(T_K+1));
        Mad = midrad(Mad_mid,Mad_rad);
      end
      % size of matrix at current iteration depth
      sizeMad = sizeMad.*sizeM(:,m);
    end
    for m=M(1)+1:M(1)+M(2)
      % we do this recursively
      % Chebyshev
      Na=sizea(m+2)-1;
      K=Nout(:,m);
      NK1=min(Na,K(1));
      NK2=min(Na,K(2));
      N2K2=min(NK1,K(1)+K(2));
      T_K = toeplitz([1:NK1+1 zeros(1,K(1)-NK1)]',[1:NK2+1 zeros(1,K(2)-NK2)]);
      % also the hankel part
      H_K = hankel([1:NK1+1 zeros(1,K(1)-NK1)]',[NK1+1:N2K2+1 zeros(1,K(2)-(N2K2-NK1))]);
      H_K(:,1)=0;
      Mad = reshape(Mad,[sizeMad(1:2)',sizea(2+m),prodsizea(m)]);
      Madzero = altzeros([size(Mad,1:2),1,size(Mad,4)],a(1));
      Mad = permute([permute(Madzero,[3 1 2 4]);permute(Mad,[3 1 2 4])],[2 3 1 4]);
      % Mad = cat(3,Madzero,Mad);
      if ~intervals
        Mad = mat2cell(Mad, sizeMad(1), sizeMad(2), ones(Na+2,1), size(Mad,4));
        Mad = cell2mat(Mad(T_K+1))+cell2mat(Mad(H_K+1));
      else
        Mad_mid = mat2cell(mid(Mad), sizeMad(1), sizeMad(2), ones(Na+2,1), size(Mad,4));
        Mad_mid_T = cell2mat(Mad_mid(T_K+1));
        Mad_mid_H = cell2mat(Mad_mid(H_K+1));
        Mad_rad = mat2cell(rad(Mad), sizeMad(1), sizeMad(2), ones(Na+2,1), size(Mad,4));
        Mad_rad_T = cell2mat(Mad_rad(T_K+1));
        Mad_rad_H = cell2mat(Mad_rad(T_K+1));
        % add rigorously
        Mad = midrad(Mad_mid_T,Mad_rad_T)+midrad(Mad_mid_H,Mad_rad_H);
      end
      sizeMad = sizeMad.*sizeM(:,m);
    end
    for m=M(1)+M(2)+1:M(1)+M(2)+M(3)
      % we do this recursively
      % Taylor
      Na=sizea(m+2)-1;
      K=Nout(:,m);
      NK1=min(Na,K(1));
      % lower triangular
      T_K = toeplitz([1:NK1+1 zeros(1,K(1)-NK1)]',[1 zeros(1,K(2))]);
      Mad = reshape(Mad,[sizeMad(1:2)',sizea(2+m),prodsizea(m)]);
      Madzero = altzeros([size(Mad,1:2),1,size(Mad,4)],a(1));
      Mad = permute([permute(Madzero,[3 1 2 4]);permute(Mad,[3 1 2 4])],[2 3 1 4]);
      % Mad = cat(3,Madzero,Mad);
      if ~intervals
        Mad = mat2cell(Mad, sizeMad(1), sizeMad(2), ones(Na+2,1), size(Mad,4));
        Mad = cell2mat(Mad(T_K+1));
      else
        Mad_mid = mat2cell(mid(Mad), sizeMad(1), sizeMad(2), ones(Na+2,1), size(Mad,4));
        Mad_mid = cell2mat(Mad_mid(T_K+1));
        Mad_rad = mat2cell(rad(Mad), sizeMad(1), sizeMad(2), ones(Na+2,1), size(Mad,4));
        Mad_rad = cell2mat(Mad_rad(T_K+1));
        Mad = midrad(Mad_mid,Mad_rad);
      end
      sizeMad = sizeMad.*sizeM(:,m);
    end
    Ta(d,dd,:,:)=Mad;
  end
end

% final reshaping
Ta=reshape(Ta,[sizea(1:2),sizeM(1,:),sizeM(2,:)]);
if nargout==2
  Ma=tensortomatrix(Ta,M);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function z = altzeros(n,intvaltest)
% Z = altzeros(N,INTVALTEST)
% produces a tensor of zeros of double or interval type
% depending on the type of INTVALTEST
% You must (really: must) bracket the dimensions in N, i.e.
% altzeros([3,4],intval('3.1'))
% The default is no intervals (Z is of type double if INTVALTEST is not provided)

if length(n)==1
  disp('please check code; possible problem in altzeros; you have to bracket dimensions')
end

z = zeros(n);
if exist('intvaltest','var') && exist('intval','file') && isintval(intvaltest(1))
  z=intval(z);
end

end
