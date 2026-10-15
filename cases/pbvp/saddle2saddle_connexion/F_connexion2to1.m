function [F,Y] = F_connexion2to1(x,params,N,L,v1u,v2u,Vs_at_ta,tau,delta,ta)

n = length(x)-1;
phi = x(1:2);
theta = x(3:n+1);

F = zeros(n+1,1);

% Locking condition in the 2D unstable manifold
F(1) = phi(1)^2 + phi(2)^2 - delta^2;

% Initial condition for the integrator
a0 = phi(1)*v1u + phi(2)*v2u;

[~,Y] = ode45(@(t,a) OK_ODE(t,a,params,N,L),[0,tau],a0);
x_at_tau = ta+Vs_at_ta*theta;
F(2:n+1) = Y(end,:)'-x_at_tau;

end
