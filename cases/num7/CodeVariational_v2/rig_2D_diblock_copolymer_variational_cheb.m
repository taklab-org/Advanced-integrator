function [Psi, r, Y_0, Z_0, Z_1] = rig_2D_diblock_copolymer_variational_cheb(ba, m, Nop_coeff, h, epsilon, sigma, nu, L, lambda)
% RIG_2D_DIBLOCK_COPOLYMER_VARIATIONAL_CHEB
%   2D diblock copolymer (Ohta-Kawasaki) variational verification
%   using upgraded Chebyshev integration method.
%
%   [Psi, r, Y_0, Z_0, Z_1] = rig_2D_diblock_copolymer_variational_cheb(
%       ba, m, Nop_coeff, h, epsilon, sigma, nu, L, lambda)
%
%   Inputs:
%     ba          : [n_cheb, m1, m2] approximate PDE solution
%     m           : [m1, m2] spatial mode sizes
%     Nop_coeff   : nonlinear operator coefficients
%     h           : time step
%     epsilon     : interface thickness parameter
%     sigma       : long-range interaction parameter
%     nu          : Chebyshev weight parameter
%     L           : [L1, L2] spatial frequency scaling
%     lambda      : (optional) additional constraint parameter
%
%   Outputs:
%     Psi   : fundamental matrix [n_cheb*m1*m2, m1*m2]
%     r     : maximum radius
%     Y_0   : Y bound per column [1, m1*m2]
%     Z_0   : Z_0 bound
%     Z_1   : Z_1 bound

if ~exist('lambda','var')
    lambda = [];
end

% Create closure for mu_k_dcCH: capture epsilon and sigma
mu_k_func = @(N, ~, L_arg) mu_k_dcCH(N, epsilon, sigma, L_arg);

% Pack epsilon and sigma as lambda parameter (not used by mu_k_func due to closure)
[Psi, r, Y_0, Z_0, Z_1] = rig_variational_cheb(ba, m, Nop_coeff, h, [epsilon, sigma], nu, L, mu_k_func);

end
