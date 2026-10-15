function [mag,r] = chebmag(a,I)% return mag (max of abs value) of Chebyshev
% a: Two-sided Chebyshev coefficients
% I: domain of the Chebyshev
%
% NOTE :: available for interval inputs, e.g., chebmag(intval(a),I)).
% For interval input, mag is a rigorous upper bound and r covers possible
% maximizers (the whole domain if coefficient uncertainty is separated).
%
arguments
    a
    I = [-1,1];
end
% dimension of Chebyshev interpolations
ndim = size(a,2);
%
if exist('intval','file') && isintval(a(1))
    mag = intval(zeros(ndim,1));
    r = intval(zeros(ndim,1));
    % |T_k|<=1 makes sum|a_k| a valid bound for every column. Columns that
    % are negligible against the largest one keep this bound: their root
    % covers stall at the subdivision limit near rounding level (PBVP data
    % near the origin) while changing the result by < 1e-8 relative.
    coefficient_bound = sup(sum(abs(a),1));
    if isreal(a)
        skip = coefficient_bound <= 1e-8*max(coefficient_bound);
    else
        skip = false(1,ndim);
    end
    domain = infsup(inf(intval(I(1))),sup(intval(I(2))));
    for i=1:ndim
        if skip(i)
            mag(i) = intval(coefficient_bound(i)); r(i) = domain;
            continue
        end
        [m1,r1] = chebmag_main(a(:,i),I);
        mag(i) = m1; r(i) = r1;
    end
else
    mag = zeros(ndim,1);
    r = zeros(ndim,1);
    for i=1:ndim
        [m1,r1] = chebmag_main(a(:,i),I);
        mag(i) = m1; r(i) = r1;
    end
end
end
%
function [maxima,rmax] = chebmag_main(a,I)
while size(a,1)>1
    if isa(a,'intval'), zero=inf(a(end))==0 && sup(a(end))==0;
    else, zero=a(end)==0;end
    if ~zero,break;end
    a(end)=[];
end
if isempty(a)
    maxima=0;rmax=I(1);return;
end
if size(a,1)==1
    if isa(a,'intval'), maxima=intval(mag(a));else,maxima=abs(a);end
    if isa(a,'intval'),rmax=infsup(inf(intval(I(1))),sup(intval(I(2))));else,rmax=I(1);end
    return;
end
if isa(a,'intval')
    % Bound coefficient uncertainty separately. Root covers then need only
    % cover the midpoint polynomial, rather than a broad family of roots.
    domain=infsup(inf(intval(I(1))),sup(intval(I(2))));
    center=intval(mid(a));eta=sum(abs(a-center));coefficient_bound=sum(abs(a));
    if ~isreal(a)
        maxima=intval(sup(coefficient_bound));rmax=domain;return;
    end
    while size(center,1)>1 && inf(center(end))==0 && sup(center(end))==0
        center(end)=[];
    end
    if size(center,1)==1
        maxima=intval(min(sup(coefficient_bound),sup(abs(center)+eta)));
        rmax=domain;return;
    end
    middle=mid((intval(I(1))+intval(I(2)))/2);
    probe=chebvals(center,intval([I(1);middle;I(2)]),I);
    lower=max(inf(abs(probe)));
    if sup(coefficient_bound)-lower<=4*sup(eta)
        maxima=intval(sup(coefficient_bound));rmax=domain;return;
    end
    roots=chebextrema(center,I);points=[roots;intval(I(:))];
    values=chebvals(center,points,I);
    upper=max(mag(values));
    maxima=intval(min(sup(coefficient_bound),sup(intval(upper)+eta)));
    rmax=domain;return;
end
% fvals at endpoints
epts = chebendpoints(a);
%
% extrema of f
r = chebextrema(a,I); rcand = [r;I'];
candidate = [chebvals(a,r,I);epts];
[maxima,indmax] = max(abs(candidate));
rmax = rcand(indmax);
end
