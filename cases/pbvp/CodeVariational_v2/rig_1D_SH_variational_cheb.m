function [Psi, r, Y_0, Z_0, Z_1] = rig_1D_SH_variational_cheb(ba, m, Nop_coeff, h, params, nu, L)
% RIG_1D_SH_VARIATIONAL_CHEB  1D Swift-Hohenberg variational verification
%                               using upgraded Chebyshev integration method.
%
%   [Psi, r, Y_0, Z_0, Z_1] = rig_1D_SH_variational_cheb(ba, m, Nop_coeff, h, params, nu, L)
%
%   Inputs:
%     ba          : [n_cheb, m] approximate PDE solution
%     m           : number of spatial modes (scalar)
%     Nop_coeff   : nonlinear operator coefficients
%     h           : time step
%     params      : [beta1, beta2] for 1D SH eigenvalues
%     nu          : Chebyshev weight parameter
%     L           : spatial frequency scaling (scalar)
%
%   Outputs:
%     Psi   : fundamental matrix [n_cheb*m, m]
%     r     : maximum radius
%     Y_0   : Y bound per column [1, m]
%     Z_0   : Z_0 bound
%     Z_1   : Z_1 bound

[Psi, r, Y_0, Z_0, Z_1] = rig_variational_cheb(ba, m, Nop_coeff, h, params, nu, L, @mu_k_1DSH);

end
