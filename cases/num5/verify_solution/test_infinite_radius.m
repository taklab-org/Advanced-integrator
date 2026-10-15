function test = test_infinite_radius(rho,W,Yforcing,Vq,anorm)
% Terminal test conditional on a VERIFIED incoming bound rho.
Y=intval(W)*intval(rho)+intval(Yforcing);
test=test_cubic_radius(Y,Vq,anorm);
end
