function [F,Y] = F_connexion4to3(x,params,N,L,v1u,v2u,v3u,v4u,Vs_at_ta,tau,delta,ta)

n = length(x)-1;
phi = x(1:4);
theta = x(5:n+1);

F = zeros(n+1,1);

% Locking condition in the 2D unstable manifold
F(1) = norm(phi,2)^2 - delta^2;

% Initial condition for the integrator
a0 = phi(1)*v1u + phi(2)*v2u + phi(3)*v3u + phi(4)*v4u;

% [~,Y] = ode45(@(t,a) OK_ODE(t,a,params,N,L),[0,tau],a0);
% x_at_tau = ta + Vs_at_ta*theta;
% F(2:n+1) = Y(end,:)'-x_at_tau;

opts = odeset('abstol',1e-16,'reltol',2.22045e-14);
% Y = ode45(@(t,a) OK_ODE(t,a,params,N,L),[0,tau],a0,opts);
Y = ode113(@(t,a) OK_ODE(t,a,params,N,L),[0,tau],a0,opts);
x_at_tau = ta+Vs_at_ta*theta;
F(2:n+1) = Y.y(:,end)-x_at_tau;

end
