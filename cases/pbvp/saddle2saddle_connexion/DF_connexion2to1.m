function DF = DF_connexion2to1(x,params,N,L,v1u,v2u,Vs_at_ta,tau,delta,ta)

m = length(x);
DF = zeros(m);
h = 1e-6; % finite difference parameter
I = eye(m);

for j = 1:m
    x_plus = x + h*I(:,j);
    x_minus = x - h*I(:,j);
    F1 = F_connexion2to1(x_plus,params,N,L,v1u,v2u,Vs_at_ta,tau,delta,ta);
    F2 = F_connexion2to1(x_minus,params,N,L,v1u,v2u,Vs_at_ta,tau,delta,ta);
    DF(:,j) = (F1-F2)/(2*h);
end

end
