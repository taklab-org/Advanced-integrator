function [x,Y] = newton_connexion2to1(x,params,N,L,v1u,v2u,Vs_at_ta,tau,delta,ta)

tol = 1e-12;

F = F_connexion2to1(x,params,N,L,v1u,v2u,Vs_at_ta,tau,delta,ta);

disp(['Before Newton ||F|| = ',num2str(norm(F,inf))])

k=0;
while (k<=40) && (norm(F)> tol)
    DF = DF_connexion2to1(x,params,N,L,v1u,v2u,Vs_at_ta,tau,delta,ta);
    x = x - DF\F;
    [F,Y] = F_connexion2to1(x,params,N,L,v1u,v2u,Vs_at_ta,tau,delta,ta);
    disp(['||F|| = ',num2str(norm(F,inf))])
    k=k+1;
end

disp(['After Newton, ||inv(DF)|| = ',num2str(norm(inv(DF)))])


end
