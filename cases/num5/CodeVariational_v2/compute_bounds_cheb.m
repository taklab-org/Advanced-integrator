function [r, Y, Z_1, Z_0, details] = compute_bounds_cheb(mu_k_sp, Nop_coeff, a_bar, h, nu, L_sp, m, build_DQ_int_func, use_shift)
% Infinite Chebyshev operator bounds with certified finite-section remainders.
%  * P K P is evaluated on a window containing the whole support of Lambda*DQ*e_j,
%    with the finite-section remainder of L^{-1} enclosed explicitly.
%  * Z_inf and Z_leN use row-side tail norms of L^{-1}Lambda.
%  * Z_A0 uses the column range that DQ can reach from the tail.
if ~exist('build_DQ_int_func','var')||isempty(build_DQ_int_func), build_DQ_int_func=@build_DQ_full_int; end
if ~exist('use_shift','var')||isempty(use_shift), use_shift=true; end
assert(use_shift,'AI:ShiftRequired','Only the shifted Ohta-Kawasaki path is certified.');
assert(nu==1,'AI:Weight','Only weight one is certified.');
assert(isequal(size(Nop_coeff),[3,4]) && Nop_coeff(3,4)==-1 && nnz(Nop_coeff)==1, ...
 'AI:Nonlinearity','Only the Ohta-Kawasaki cubic is supported.');
assert(all(rad(intval(a_bar))==0,'all'),'AI:Center','An exact polynomial center is required.');
h_iv=intval(h); L_iv=intval(L_sp);
Ma=size(a_bar); n=Ma(end); N2=n-1; num_sp=prod(m); dimP=num_sp*n; Mc=[m,n];
lam=intval(mu_k_sp(:)); lam_sup=intval(sup(max(lam)));
mu=h_iv/2*(lam_sup-lam); mumax=sup(max(abs(mu)));
A=reshape(abs(mid(a_bar)),[],n); na=find(max(A,[],1)>0,1,'last')-1; deg=2*na;
assert(N2>=deg,'Chebyshev window too small for the convolution support.');
Md=N2+deg+2;                               % DQ*e_j has support in 0..Md
Mw=max(Md,ceil(8*mumax)+40);               % window for L^{-1}: indices 0..Mw
aw=zeros([Ma(1:end-1),Md+1]); idx=repmat({':'},1,numel(Ma)); idx{end}=1:n; aw(idx{:})=mid(a_bar);
DQw=build_DQ_int_func(intval(reshape(mu_k_sp,m)),intval(Nop_coeff),intval(aw),h_iv,Md+1,lam_sup,L_iv);
w=2*intval(ones(Mw+1,1)); w(1)=intval(1);
ws=spatial_weights(m); om=kron(ws,[intval(1);2*intval(ones(n-1,1))]);
colsP=reshape((0:num_sp-1)*(Md+1)+(1:n).',1,[]);
Dcols=DQw(:,colsP); clear DQw; omw=kron(ws,w(1:Md+1));
D=intval(0); for j=1:dimP, D=max(D,sum(omw.*abs(Dcols(:,j)))/om(j)); end
D=intval(sup(D));
Lam=intval(conversion_op(Mw));
KPP=intval(zeros(dimP)); dcol=intval(zeros(1,dimP));
G=cell(num_sp,1); rho=intval(zeros(num_sp,1)); rho1=rho;
for s=1:num_sp
  [rho(s),G{s}]=compute_rho_k(Mw,mu(s));
  assert(sup(rho(s))<1,'rho >= 1 on the window.');
  g0=sum(w.*abs(G{s}(:,1))); rho1(s)=(g0+1)/intval(Mw+1);
  GL=G{s}*Lam;
  X=GL(:,1:Md+1)*Dcols((s-1)*(Md+1)+(1:Md+1),:);
  KPP((s-1)*n+(1:n),:)=h_iv/2*X(1:n,:);
  dcol=dcol+ws(s)*(abs(mu(s))*rho1(s)/(1-rho(s)))*abs(X(end,:));
end
DF=KPP; for i=1:dimP, DF(i,i)=DF(i,i)+1; end
A0=inv(mid(DF)); A0i=intval(A0);
B=-A0i*DF; for i=1:dimP, B(i,i)=B(i,i)+1; end
A0n=wnorm_op(A0i,om);
Z_0=sup(wnorm_op(B,om)+A0n*(h_iv/2)*max(dcol./om.'));
% tail norms of L^{-1}Lambda, mode by mode on a window adapted to mu
T=intval(0); Tp=intval(0);
for s=1:num_sp
  % Bound each interval eigenvalue separately, without midpoint deduplication.
  [ts,tp]=tail_norms(mu(s),N2,deg,mumax); T=max(T,ts); Tp=max(Tp,tp);
end
Z_inf=h_iv/2*T*D; Z_leN=Z_inf; Z_A0=A0n*(h_iv/2)*Tp*D;
t1=Z_inf^2+(1+Z_inf)*Z_leN*Z_A0+Z_A0*Z_inf^2;
t2=intval(Z_0)*(1+Z_leN*(1+Z_inf))+Z_leN*Z_A0*(1+Z_inf);
Z_1=sup(max(t1,t2));
details=struct('Z_0',Z_0,'Z_1',Z_1,'Z_inf',sup(Z_inf),'Z_leN',sup(Z_leN),'Z_A0',sup(Z_A0), ...
 'rho_vals',sup(rho),'mu_paper',mu,'A0',A0,'A0n',A0n,'Mc',Mc,'Nop_coeff',intval(Nop_coeff), ...
 'a_bar',intval(a_bar),'L_sp',L_iv,'lambda_sup',lam_sup,'build_DQ_int_func',build_DQ_int_func, ...
 'n_fourier',Ma(1:end-1),'Mw',Mw,'Md',Md,'Dcols',Dcols,'G',{G},'rho',rho,'rho1',rho1,'Lam',Lam);
Y=nan; r=nan;

end
function v=wnorm_op(B,om)
v=intval(0); for j=1:size(B,2), v=max(v,sum(om.*abs(B(:,j)))/om(j)); end
end
function ws=spatial_weights(m)
[k1,k2]=ndgrid(0:m(1)-1,0:m(2)-1); ws=intval(2.^((k1(:)>0)+(k2(:)>0)));
end
function [T,Tp]=tail_norms(mu,N2,deg,mumax)
Mt=max(N2+40,ceil(8*mumax)+40);
[rho,G]=compute_rho_k(Mt,mu); assert(sup(rho)<1);
w=2*intval(ones(Mt+1,1)); w(1)=intval(1);
g0=sum(w.*abs(G(:,1))); r1=(g0+1)/intval(Mt+1);
R=G*intval(conversion_op(Mt));
T=intval(0); Tp=intval(0); j0=max(0,N2+1-deg);
for j=0:Mt-1
  x=R(:,j+1); dl=abs(mu)*r1*abs(x(end))/(1-rho);
  T=max(T,(sum(w(N2+2:end).*abs(x(N2+2:end)))+dl)/w(j+1));
  if j>=j0, Tp=max(Tp,(sum(w(1:N2+1).*abs(x(1:N2+1)))+dl)/w(j+1)); end
end
% columns Mt, Mt+1 and beyond: explicit A_k*Lambda*e_j, then the Neumann factor
cM =(sum(w.*abs(G(:,Mt)))+(1+g0)/intval(Mt+1))/2;
cM1=(sum(w.*abs(G(:,Mt+1)))+(1+g0)/intval(Mt+2))/2;
cF =(1/intval(Mt+1)+1/intval(Mt+3)+2*g0/(intval(Mt+1)*intval(Mt+3)))/2;
c=max(max(cM,cM1),cF)/(1-rho);
T=max(T,c); Tp=max(Tp,c);
end
