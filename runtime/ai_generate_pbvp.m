function ai_generate_pbvp(out)
% Rebuild the numerical candidate, not a certificate, from a mode seed.
% No saved flow, eigenbasis, inverse, or interval bounds are imported.
params=[0;1;.39];N=[12,12];L=[1;1.1];tau=.1;delta=.01;udim=4;
n=prod(N)-1;
% Origin linearization is diagonal. Fix the four coordinate eigenvectors
% explicitly to avoid eig ordering/sign differences across LAPACK versions.
V=zeros(n,4);idx=[2,14,13,24];sgn=[-1,1,-1,1];
for k=1:4,V(idx(k),k)=sgn(k);end
v1u=V(:,1);v2u=V(:,2);v3u=V(:,3);v4u=V(:,4);
A=zeros(N);A(3,2)=.1;ta=A(:);ta(1)=[];
ta=newton(ta,params,N,L);assert(norm(F_DC_2D(ta,params,N,L),inf)<1e-11);
[Q,D]=eig(DF_DC_2D(ta,params,N,L));ev=diag(D);
assert(all(abs(imag(ev))<1e-8)&&sum(real(ev)>0)==3);
Vs_at_ta=real(Q(:,real(ev)<0));Vs_at_ta=fliplr(Vs_at_ta);
assert(size(Vs_at_ta,2)==n-3);
theta=[0;delta;0;0];opts=odeset('AbsTol',1e-16,'RelTol',2.22045e-14);
Y=ode113(@(t,a)OK_ODE(t,a,params,N,L),[0,tau],V*theta,opts);
phi=Vs_at_ta\(Y.y(:,end)-ta);x=[theta;phi];
% Refine all boundary coordinates with a fresh Newton solve. A full step
% is accepted only when it decreases the actual boundary residual.
for iter=1:12
 [F,Y]=F_connexion4to3(x,params,N,L,v1u,v2u,v3u,v4u,Vs_at_ta,tau,delta,ta);
 if norm(F,inf)<2e-14,break;end
 J=DF_connexion4to3(x,params,N,L,v1u,v2u,v3u,v4u,Vs_at_ta,tau,delta,ta);
 dx=J\F;alpha=1;
 for ls=1:20
  candidate=x-alpha*dx;Fc=F_connexion4to3(candidate,params,N,L,v1u,v2u,v3u,v4u,Vs_at_ta,tau,delta,ta);
  if norm(Fc,inf)<norm(F,inf),break;end
  alpha=alpha/2;
 end
 x=candidate;
end
[F,Y]=F_connexion4to3(x,params,N,L,v1u,v2u,v3u,v4u,Vs_at_ta,tau,delta,ta);
assert(norm(F,inf)<1e-11,'Numerical PBVP candidate failed; no proof claimed.');
save(fullfile(out,'pathdata4to3.mat'),'params','N','L','tau','delta','udim', ...
 'v1u','v2u','v3u','v4u','Vs_at_ta','ta','x','Y','-v7.3');
end
