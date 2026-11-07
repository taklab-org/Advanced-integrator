function [c, k] = newton_cheb(L_full, Lam_full, DQ_full, beta, h, c0)
% NEWTON_CHEB  Newton iteration for the Chebyshev integration formulation.
%
%   [c, k] = newton_cheb(L_full, Lam_full, DQ_full, beta, h, c0)
%
%   Solves: G(c) = L*c + (h/2)*Lambda*DQ*c - beta = 0
%   Since G is linear in c, Newton converges in one step.
%
%   Inputs:
%     L_full, Lam_full, DQ_full : operators (matrices)
%     beta     : right-hand side
%     h        : time step
%     c0       : initial guess
%
%   Outputs:
%     c : solution vector
%     k : number of iterations

tol = 5e-11;
total_dim = length(c0);
c = c0;

% Build the constant Jacobian
DF = DF_cheb(L_full, Lam_full, DQ_full, h);

% Evaluate F
F = F_cheb(c, L_full, Lam_full, DQ_full, beta, h);

k = 0;
k_max = 10;  % should converge in 1 step for linear problem

while k <= k_max && norm(F, 1) > tol
    c = c - DF \ F;
    F = F_cheb(c, L_full, Lam_full, DQ_full, beta, h);
    k = k + 1;
end

if norm(F, 1) > tol
    error('Newton did not converge: ||F|| = %e after %d iterations', norm(F,1), k);
end

disp(['||F(c)|| = ', num2str(norm(F,1)), ', # Newton iterations = ', num2str(k)]);

end
