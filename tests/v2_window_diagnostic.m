function R=v2_window_diagnostic(ba,h,m,params,L,n_pad,extra,tag,doY)
% Audit harness for CodeVariational_v2: code bounds against large-window operator norms.
% Double-precision truth values are diagnostics, not proofs.
if nargin<9, doY=true; end
ep=params(3); sigma=params(2);
ba=mid(intval(ba)); [n0,N1,N2f]=size(ba); n=n0+n_pad; N2=n-1;
e=zeros(n,N1,N2f); e(1:n0,:,:)=ba; a_bar=permute(e,[2,3,1]);
Nop=[0,0,0,0;0,0,0,0;0,0,0,-1]; num_sp=prod(m);
cd0=zeros([m,n]); muf=mu_k_dcCH(size(cd0)-1,ep,sigma,L); mu_sp=reshape(muf(:,:,1),[],1);
t0=tic;
evalc('[~,~,Z1c,Z0c,det]=compute_bounds_cheb(mu_sp,Nop,a_bar,h,1,L,m,@build_DQ_full_int,true);');
tcode=toc(t0);
lam=max(mu_sp); mup=(h/2)*(lam-mu_sp);
[~,b]=powerconvcos(ba,2); deg=size(b,1)-1;
M=n+deg+extra; dim=num_sp*M;
DQ=assembleDQ(b,m,M,L,lam);
Lam=conversion_op(M-1); om=omega(m,M);
Kb=zeros(dim); LLq=zeros(num_sp,1); Tn=zeros(num_sp,1);
nidx=repmat((0:M-1).',num_sp,1);
for s=1:num_sp
  id=(s-1)*M+(1:M);
  G=complete_op(M-1,mup(s))\Lam;
  Kb(id,:)=(h/2)*(G*DQ(id,:));
  w=[1;2*ones(M-1,1)];
  cols=1:(M-5);
  tailrows=(N2+2):M;
  Tn(s)=max((w(tailrows).'*abs(G(tailrows,cols)))./w(cols).');
end
P=find(nidx<=N2); Q=find(nidx>N2);
vq=nidx(Q)<=M-deg-6;       % columns not cut by the window
Kpp=Kb(P,P); Kpq=Kb(P,Q); Kqp=Kb(Q,P); Kqq=Kb(Q,Q);
A0=det.A0; wP=om(P); wQ=om(Q);
% reorder: code uses block order on n, so P inherits it from the window
nI=numel(P);
Bp=eye(nI)-A0*(eye(nI)+Kpp);
tZ0=opn(Bp,wP,wP);
tZinf=opn(Kqq(:,vq),wQ,wQ(vq));
tZle=opn(Kqp,wQ,wP);
tZA=opn(A0*Kpq(:,vq),wP,wQ(vq));
AK=A0*Kpq; QK=eye(numel(Q))-Kqq; KqpB=Kqp*Bp;
top=Bp+AK*(QK*Kqp); bot=-(QK*KqpB);
cP=max((wP.'*abs(top)+wQ.'*abs(bot))./wP.');
K2=Kqq*Kqq;
top2=-AK*K2; bot2=K2+QK*(Kqp*AK);
cQ=max((wP.'*abs(top2(:,vq))+wQ.'*abs(bot2(:,vq)))./wQ(vq).');
tZ1=max(cP,cQ);
tD=opn(DQ(:,nidx<=M-deg-6),om,om(nidx<=M-deg-6));
R=struct('tag',tag,'n',n,'n0',n0,'deg',deg,'h',h,'mumax',max(mup),'rho',max(det.rho_vals), ...
 'Z0',det.Z_0,'tZ0',tZ0,'Zinf',det.Z_inf,'tZinf',tZinf,'Zle',det.Z_leN,'tZle',tZle, ...
 'ZA',det.Z_A0,'tZA',tZA,'Z1',det.Z_1,'tZ1',tZ1,'D',NaN,'tD',tD,'T',max(Tn),'tcode',tcode);
fprintf('%-20s n=%3d n0=%3d deg=%3d h=%.5f mu_max=%.3g rho=%.3g  (code %.0fs)\n',tag,n,n0,deg,h,R.mumax,R.rho,tcode);
fprintf('   quantity    code bound    window truth   code/truth\n');
pr('Z_0',det.Z_0,tZ0); pr('Z_inf',det.Z_inf,tZinf); pr('Z_<=N',det.Z_leN,tZle); pr('Z_A0',det.Z_A0,tZA); pr('Z_1',det.Z_1,tZ1);
if doY
  [Lf,Lamf,~]=build_Lk_full(mu_sp,h,n,true);
  mu_md=reshape(mu_sp,m);
  evalc('DQf=build_DQ_full(mu_md,Nop,a_bar,h,n,lam,L);');
  DFf=Lf+(h/2)*Lamf*DQf;
  Lb=zeros(dim); Lmb=zeros(dim);
  for s=1:num_sp, id=(s-1)*M+(1:M); Lb(id,id)=complete_op(M-1,mup(s)); Lmb(id,id)=Lam; end
  Dbig=Lb+(h/2)*Lmb*DQ;
  Bt=zeros(dim,num_sp); for s=1:num_sp, Bt((s-1)*M+1,s)=1; end
  Ut=Dbig\Bt;
  Yc=zeros(1,num_sp); Y0c=Yc; Y3c=Yc; err=Yc; rc=Yc;
  den=1-det.Z_0-det.Z_1;
  for i=1:num_sp
    beta=zeros(num_sp*n,1); beta((i-1)*n+1)=1;
    c=DFf\beta;
    [Yi,Y0i,Y3i]=compute_Y_cheb(c,beta,det,1,h);
    Yc(i)=Yi; Y0c(i)=Y0i; Y3c(i)=Y3i; rc(i)=Yi/den;
    Ub=zeros(dim,1); Ub(P)=c;
    err(i)=sum(om.*abs(Ut(:,i)-Ub));
  end
  R.r=rc; R.err=err; R.Y=Yc; R.Y0=Y0c; R.Y3=Y3c;
  [worst,iw]=max(err./rc);
  fprintf('   columns: max r1=%.3e  max err=%.3e  worst err/r1=%.3g (col %d: err=%.3e r1=%.3e, Y0=%.3e Y3=%.3e)\n', ...
    max(rc),max(err),worst,iw,err(iw),rc(iw),Y0c(iw),Y3c(iw));
end
end
function pr(s,c,t)
fprintf('   %-9s %.5e   %.5e   %8.3g %s\n',s,c,t,c/t,repmat('<-- below truth',1,c<t));
end
function v=opn(B,wr,wc)
v=max((wr.'*abs(B))./wc.');
end
function w=omega(m,M)
[k1,k2]=ndgrid(0:m(1)-1,0:m(2)-1); ws=2.^((k1(:)>0)+(k2(:)>0));
wn=[1;2*ones(M-1,1)];
w=kron(ws,wn);
end
function DQ=assembleDQ(b,m,M,L,lam)
num_sp=prod(m);
[k1,k2]=ndgrid(0:m(1)-1,0:m(2)-1); k=[k1(:),k2(:)];
K2=(k(:,1)*L(1)).^2+(k(:,2)*L(2)).^2;
[Lc,S1,S2]=size(b);
[E1,E2]=ndgrid(0:M-1,0:M-1); I1=abs(E1-E2); I2=abs(E1+E2); ok1=I1<Lc; ok2=I2<Lc;
DQ=zeros(num_sp*M);
for s1=1:num_sp
 if K2(s1)==0, continue; end
 for s2=1:num_sp
  blk=zeros(M);
  if k(s2,1)>0, g1=[1,-1]; else, g1=1; end
  if k(s2,2)>0, g2=[1,-1]; else, g2=1; end
  for a=g1
   for c=g2
    p=abs(k(s1,1)-a*k(s2,1)); q=abs(k(s1,2)-c*k(s2,2));
    if p<S1 && q<S2
     col=b(:,p+1,q+1);
     v=zeros(M); v(ok1)=col(I1(ok1)+1); blk=blk+v;
     v=zeros(M); v(ok2)=col(I2(ok2)+1); v(:,1)=0; blk=blk+v;
    end
   end
  end
  DQ((s1-1)*M+(1:M),(s2-1)*M+(1:M))=-3*K2(s1)*blk;
 end
end
DQ=DQ+lam*eye(num_sp*M);
end
