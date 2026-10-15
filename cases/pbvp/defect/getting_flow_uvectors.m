function UVu = getting_flow_uvectors(params,N,L,Y,v1u,v2u,tau)
%
opts = odeset('abstol',1e-16,'reltol',2.22045e-14);
Uv1u = ode113(@(t,a) OK_VarODE(t,a,params,N,L,Y),[0,tau],v1u,opts);
Uv2u = ode113(@(t,a) OK_VarODE(t,a,params,N,L,Y),[0,tau],v2u,opts);
%
UVu = [Uv1u.y(:,end), Uv2u.y(:,end)];
%
end
