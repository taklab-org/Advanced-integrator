function cout = chebdiff(c)% return differentiated Two-sided Chebyshev
% c: Two-sided Chebyshev coefficients
%
% NOTE :: available for interval inputs
%
  [n, m] = size(c);

  if n<=1
    if isa(c,'intval'), cout=intval(zeros(0,m)); else, cout=zeros(0,m); end
    return;
  end

  if exist('intval','file') && isintval(c(1))
    cout = intval(zeros(n-1, m));
    w = intval(repmat(2*(1:n-1)', 1, m));
  else
    cout = (zeros(n-1, m));
    w = repmat(2*(1:n-1)', 1, m);
  end

  v = w.*c(2:end,:);
  cout(n-1:-2:1,:) = vcumsum(v(n-1:-2:1,:));
  cout(n-2:-2:1,:) = vcumsum(v(n-2:-2:1,:));
  cout(1,:) = .5*cout(1,:);
end

function csum = vcumsum(V)
  n=size(V,1);
  csum=tril(ones(n))*V;
end
