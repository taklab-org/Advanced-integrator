function local = getting_bridge_bounds(seg,h,params,L,nu,frame)
% Full-PDE one-step bounds in a fixed real finite coordinate frame.
m=frame.m;n=prod(m);D=frame.D;Di=frame.Di;V=frame.V;Vi=frame.Vi;
assert(inf(intval(nu))>=1);
[semi,semiq,b2,~,~,~,finite]=getting_int_semi_matrix(seg.ba,h,m,params,L,nu,true);
[psi,cout]=getting_coupling_vectors(b2,m,[size(seg.ba,2),size(seg.ba,3)],nu);
[v0,vq,mu]=getting_Vinf(h,[size(seg.ba,2),size(seg.ba,3)],m,params,L);
cin=abs(Vi)*(D*(abs(semiq)*(3*psi)));
co=((Di*(vq*3*cout)).')*abs(V);g=vq*3*wnorm(b2,intval(nu));
s=intval(1)-g-co*cin;
if ~(inf(s)>0),error('getting_bridge_bounds:Schur','Fixed-frame Schur condition failed.');end
B=intval(zeros(n+1));B(1:n,1:n)=intval(eye(n))+cin*(co/s);
B(1:n,end)=cin/s;B(end,1:n)=co/s;B(end,end)=1/s;
G=intval(zeros(n+1));G(1:n,end)=cin;G(end,1:n)=co;G(end,end)=g;
[supz,~]=getting_finite_flow(finite,h,Vi*D,Di*V);
[~,Phi]=getting_finite_flow(finite,h,D,Di);
base=intval(zeros(n+1));base(1:n,1:n)=supz;base(end,end)=1;
S=B*base;R=G*S;sigma=exp(mu*intval(h));
Phi_z=Vi*Phi*V;endpoint=intval(zeros(n+1));
endpoint(1:n,1:n)=abs(Phi_z);endpoint(end,end)=sigma;
H2=getting_block_norms(endpoint+R,n); % Combine columns BEFORE taking norms.
Vfin=abs(Vi)*(D*semi*Di);Vfinq=abs(Vi)*(D*semiq*Di);
qbase=intval(zeros(n+1));qbase(1:n,1:n)=Vfinq;qbase(end,end)=vq;
Vq=B*qbase;
forcing=getting_time_forcing(seg.defect,m,h,params,L,nu,B,Vfin);
local=struct('S',S,'S2',getting_block_norms(S,n),'G',G,'B',B,'R',R, ...
    'H2',H2,'Phi',Phi,'sigma',sigma,'finite',finite,'Vq',Vq, ...
    'Vq2',getting_block_norms(Vq,n),'Vq_norm',sup(norm(Vq,1)), ...
    'forcing',forcing,'schur',inf(s),'neumann',sup(norm(finite.K,inf)),'h',h);
end
