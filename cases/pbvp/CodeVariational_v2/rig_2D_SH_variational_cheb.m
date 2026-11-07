function [Psi, r, Y_0, Z_0, Z_1] = rig_2D_SH_variational_cheb(ba, m, Nop_coeff, h, lambda, nu, L)
% RIG_2D_SH_VARIATIONAL_CHEB  2D Swift-Hohenberg variational verification
%                               using upgraded Chebyshev integration method.
%
%   [Psi, r, Y_0, Z_0, Z_1] = rig_2D_SH_variational_cheb(ba, m, Nop_coeff, h, lambda, nu, L)
%
%   Inputs:
%     ba          : [n_cheb, m1, m2] approximate PDE solution
%     m           : [m1, m2] spatial mode sizes
%     Nop_coeff   : nonlinear operator coefficients
%     h           : time step
%     lambda      : physical parameter for SH eigenvalues
%     nu          : Chebyshev weight parameter
%     L           : [L1, L2] spatial frequency scaling
%
%   Outputs:
%     Psi   : fundamental matrix [n_cheb*m1*m2, m1*m2]
%     r     : maximum radius
%     Y_0   : Y bound per column [1, m1*m2]
%     Z_0   : Z_0 bound
%     Z_1   : Z_1 bound

[Psi, r, Y_0, Z_0, Z_1] = rig_variational_cheb(ba, m, Nop_coeff, h, lambda, nu, L, @mu_k_SH);

end
