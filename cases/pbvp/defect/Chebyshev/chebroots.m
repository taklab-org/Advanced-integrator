function r = chebroots(c,I)
% Roots of c_0 + c_1 T_1 + ... on I (two-sided coefficients).
% Interval input returns a cover of ALL possible real roots. Unresolved
% boxes are retained instead of being discarded by a numerical eig.
arguments
    c
    I = [-1,1];
end
isiv = isa(c,'intval') || isa(I,'intval');
if isiv, c=intval(c);end
c = c(:);
while size(c,1)>1 && exact_zero(c(end),isiv)
    c(end)=[];
end
if isiv
    assert(prod(size(I))==2 && inf(intval(I(2))-intval(I(1)))>0);
    if isempty(c) || (size(c,1)==1 && contains_zero(c(1)))
        z=infsup(-1,1); % The zero polynomial can vanish everywhere.
    elseif size(c,1)==1
        z=intval(zeros(0,1));
    elseif size(c,1)==2 && ~contains_zero(c(2))
        z=intersect_unit(-c(1)/c(2));
    else
        z=interval_root_cover(c);
    end
    r=(intval(1)-z)*intval(I(1))/2+(intval(1)+z)*intval(I(2))/2;
    return;
end
if isempty(c) || all(c==0)
    error('chebroots:ZeroPolynomial','The zero polynomial has no finite root list.');
end
n=size(c,1)-1;
if n==0
    r=zeros(0,1);
elseif n==1
    r=-c(1)/c(2);
else
    A=zeros(n);A(2,1)=1;
    for j=2:n-1
        A(j-1,j)=0.5;A(j+1,j)=0.5;
    end
    A(:,n)=-0.5*c(1:n)/c(end);
    A(n-1,n)=A(n-1,n)+0.5;
    r=eig(A);
end
htol=100*eps;
r=real(r(abs(imag(r))<htol));
r=sort(r(abs(r)<=1+htol));r=max(-1,min(1,r));
r=(1-r)*I(1)/2+(1+r)*I(2)/2;
end

function yes=exact_zero(x,isiv)
if isiv,yes=inf(x)==0 && sup(x)==0;else,yes=x==0;end
end

function yes=contains_zero(x)
yes=inf(x)<=0 && sup(x)>=0;
end

function z=intersect_unit(x)
lo=max(-1,inf(x));hi=min(1,sup(x));
if lo>hi,z=intval(zeros(0,1));else,z=infsup(lo,hi);end
end

function roots=interval_root_cover(c)
if all(inf(c)<=0 & sup(c)>=0)
    roots=infsup(-1,1);return; % The coefficient box contains the zero polynomial.
end
dc=chebdiff(c);
k=intval((0:size(c,1)-1)');
lip=sum(abs(c).*k.^2); % |T_k'| <= k^2 throughout [-1,1].
pending=[-1,1];boxes=zeros(0,2);processed=0;
maxboxes=512;tolerance=1e-12;
while ~isempty(pending)
    lo=pending(end,1);hi=pending(end,2);pending(end,:)=[];
    J=infsup(lo,hi);m=mid(J);radius=sup(abs(J-intval(m)));
    fmid=chebvals(c,intval(m));fbox=chebvals(c,J);
    meanbox=fmid+infsup(-radius,radius)*lip;
    if ~contains_zero(fbox) || ~contains_zero(meanbox),continue;end
    processed=processed+1;
    derivative=chebvals(dc,J);
    if ~contains_zero(derivative)
        newton=intval(m)-fmid/derivative;
        newlo=max(lo,inf(newton));newhi=min(hi,sup(newton));
        if newlo>newhi,continue;end
        lo=newlo;hi=newhi;J=infsup(lo,hi);
    end
    if hi-lo<=tolerance || processed>=maxboxes
        boxes(end+1,:)=[lo,hi]; %#ok<AGROW>
        if processed>=maxboxes
            boxes=[boxes;pending]; %#ok<AGROW>
            break;
        end
    else
        m=mid(J);
        if m<=lo || m>=hi
            boxes(end+1,:)=[lo,hi]; %#ok<AGROW>
        else
            pending=[pending;lo,m;m,hi]; %#ok<AGROW>
        end
    end
end
if isempty(boxes),roots=intval(zeros(0,1));return;end
boxes=sortrows(boxes,1);merged=boxes(1,:);
for j=2:size(boxes,1)
    if boxes(j,1)<=merged(end,2)
        merged(end,2)=max(merged(end,2),boxes(j,2));
    else
        merged(end+1,:)=boxes(j,:); %#ok<AGROW>
    end
end
roots=infsup(merged(:,1),merged(:,2));
end
