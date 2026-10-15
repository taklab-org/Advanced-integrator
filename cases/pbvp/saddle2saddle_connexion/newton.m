function a = newton(a,params,N,L)

tol = 1e-12;

F = F_DC_2D(a,params,N,L);

display(['Before Newton ||F|| = ',num2str(norm(F,1))])

k=0;
while (k<=40) && (norm(F)> tol)
    DF = DF_DC_2D(a,params,N,L);
    a = a - DF\F;
    F = F_DC_2D(a,params,N,L);
    display(['||F|| = ',num2str(norm(F,1))])
    k=k+1;
end

return
