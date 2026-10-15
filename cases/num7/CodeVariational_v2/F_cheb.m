function F = F_cheb(c, L_full, Lam_full, DQ_full, beta, h)
% F_CHEB  Evaluate the zero-finding map F in the paper's formulation.
%
%   F = F_cheb(c, L_full, Lam_full, DQ_full, beta, h)
%
%   The equation is: L*c + (h/2)*Lambda*DQ*c - beta = 0
%   where beta encodes the initial condition.
%
%   Equivalently in the paper's F formulation:
%   F(c) = c + (h/2)*L^{-1}*Lambda*DQ*c - L^{-1}*beta
%
%   But for Newton iteration, it is more efficient to solve:
%   G(c) = L*c + (h/2)*Lambda*DQ*c - beta = 0
%   since we don't need to invert L explicitly.
%
%   Inputs:
%     c        : coefficient vector (total_dim x 1)
%     L_full   : block-diagonal L operator
%     Lam_full : block-diagonal Lambda operator
%     DQ_full  : DQ operator matrix
%     beta     : right-hand side (initial condition)
%     h        : time step

F = L_full * c + (h/2) * Lam_full * DQ_full * c - beta;

end
