function Cinf = compute_Cinf(params,m,L,gam,xi)
sigma   = params(2);
epsilon = params(3);
%
m1 = m(1); m2 = m(2);
N1 = 10; N2 = 15;
L1 = L(1); L2 = L(2);
%
[k1,k2] = ndgrid(0:N1-1,0:N2-1);
%
bkL2 = (k1*L1).^2 + (k2*L2).^2;
mu_k = bkL2 .* (-epsilon^2*bkL2 + 1) - sigma;
%
ctemp = bkL2 ./ abs(mu_k).^gam;
ctemp(1:m1,1:m2) = 0;
% sup(ctemp)
Cinf = gam/(intval('e')*xi) * max(sup(ctemp),[],'all');
