function [Y_lin1, Y_lin2, Y_lin3, Y_lin4] = getting_flow_fromDP(params,N,L,Y,DP,tau)
%
v1u = DP(:,1);
v2u = DP(:,2);
v3u = DP(:,3);
v4u = DP(:,4);
%
opts = odeset('abstol',1e-16,'reltol',2.22045e-14);
Y_lin1 = ode113(@(t,a) OK_VarODE(t,a,params,N,L,Y),[0,tau],v1u,opts);
Y_lin2 = ode113(@(t,a) OK_VarODE(t,a,params,N,L,Y),[0,tau],v2u,opts);
Y_lin3 = ode113(@(t,a) OK_VarODE(t,a,params,N,L,Y),[0,tau],v3u,opts);
Y_lin4 = ode113(@(t,a) OK_VarODE(t,a,params,N,L,Y),[0,tau],v4u,opts);
%
% UVu = Y_lin(:,end);
% UVu = [Y_lin1.y(:,end), Y_lin2.y(:,end)];
%
end
