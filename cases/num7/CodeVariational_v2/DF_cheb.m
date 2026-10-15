function DF = DF_cheb(L_full, Lam_full, DQ_full, h)
% DF_CHEB  Compute the Jacobian of the zero-finding map.
%
%   DF = DF_cheb(L_full, Lam_full, DQ_full, h)
%
%   For the equation G(c) = L*c + (h/2)*Lambda*DQ*c - beta = 0:
%   DG = L + (h/2)*Lambda*DQ
%
%   This is constant (independent of c) since the problem is linear.

DF = L_full + (h/2) * Lam_full * DQ_full;

end
