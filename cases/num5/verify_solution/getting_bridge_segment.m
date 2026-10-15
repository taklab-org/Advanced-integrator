function seg = getting_bridge_segment(a0,N,h,m,params,L,nu)
% Numerical polynomial candidate followed by FULL interval defect evaluation.
% The ODE solver and Fourier/Chebyshev truncations provide no proof by themselves.
assert(params(1)==0 && all(N>=m));a0=a0(:);
opts=odeset('RelTol',2.5e-13,'AbsTol',1e-15, ...
    'Jacobian',@(t,y)DF_OK_2D(y,mid(params),N,mid(L)));
sol=ode15s(@(t,y)F_OK_2D(y,mid(params),N,mid(L)),[0,h],a0(2:end),opts);
for degree=[8,16,32,64]
    nodes=cos(pi*(0:degree)'/degree);times=(h/2)*(1+nodes);
    values=deval(sol,times.').';extended=[values;values(degree:-1:2,:)];
    a=real(fft(extended,[],1))/degree;a=a(1:degree+1,:);
    a([1,end],:)=a([1,end],:)/2;
    if max(abs(a(end-3:end,:)),[],'all')<=1e-13*max(1,max(abs(a),[],'all')),break;end
end
two=zeros(degree+1,prod(N));two(:,2:end)=a;
a_end=sum(two,1).';one=intval(two);one(2:end,:)=one(2:end,:)/2;
ba=reshape(one,degree+1,N(1),N(2));
du=chebdiff(intval(two))/(intval(h)/2);du(2:end,:)=du(2:end,:)/2;
du_ext=intval(zeros(3*degree+1,3*N(1)-2,3*N(2)-2));
du_ext(1:degree,1:N(1),1:N(2))=reshape(du,degree,N(1),N(2));
defect=du_ext-F_OK_2D_ext(one,intval(params),degree+1,N,intval(L));
start_error=compute_eps2(ba,intval(a0));
end_value=reshape(ba(1,:,:)+2*sum(ba(2:end,:,:),1),N);
end_error=end_value-reshape(intval(a_end),N);
seg=struct('ba',ba,'a_end',a_end,'defect',defect,'start_error',start_error, ...
    'end_error',end_error,'banorm',wnorm(reshape(chebmag_oneside(ba),N),intval(nu)));
end
