function [Y, Y_0, Y_3] = compute_Y_cheb(c, e_j, details, nu, h)
% COMPUTE_Y_CHEB  Compute the Y bound (F-formulation).
%
%   [Y, Y_0, Y_3] = compute_Y_cheb(c, e_j, details, nu, h)
%
%   F-formulation:
%     F(c) = c + (h/2)*L^{-1}*Lambda*DQ*c - L^{-1}*beta   [eq (13)]
%     Y_0 = ||A_0 * pi^{<=N2} F(c)||_{1,omega}
%     Y_3 = ||pi^{>N2} L^{-1}*(rhs)||_{1,omega}   (Chebyshev tail)
%
%   rhs = (h/2)*Lambda*DQ*c - beta has FINITE support (up to P*N2+1)
%   because U_bar and c have finite Chebyshev support.
%   Y_3 includes an extended-window evaluation and a certified
%   infinite finite-section remainder.
%
%   For finite Fourier, linear problem (Y_1=Y_2=0):
%     Y = Y_0*(1 + Z_{<=N}*(1+Z_inf)) + Y_3*(1+Z_inf)*(Z_{A_0}+1)   [eq (32)]

h_iv = intval(h);
c_iv = intval(c);
e_j_iv = intval(e_j);
Mc = details.Mc; n_cheb = Mc(end); Mc_sp = Mc(1:end-1);
num_sp = prod(Mc_sp); Mw = details.Mw; Md = details.Md;
% Certified finite-section evaluation.  F(c) = L^{-1} R with R = L c + (h/2) Lambda DQ c - beta.
% R is formed first, so the radius of G only meets the small residual.
% c has finite support, so DQ c is one matrix product (no interval FFT per column).
DC = reshape(details.Dcols*c_iv, Md+1, num_sp);
LamM = details.Lam(:,1:Md+1);
w_sp = build_spatial_weights(Mc_sp);
w = 2*intval(ones(Mw+1,1)); w(1) = intval(1);
cB = reshape(c_iv, n_cheb, num_sp);
eB = reshape(e_j_iv, n_cheb, num_sp);
PF = intval(zeros(n_cheb, num_sp));
Y_3 = intval(0); dtot = intval(0);
for s = 1:num_sp
    mu_s = details.mu_paper(s);
    Rs = (h_iv/2)*(LamM*DC(:,s));
    Rs(1:n_cheb) = Rs(1:n_cheb) - eB(:,s);
    Rs(1) = Rs(1) + cB(1,s) + sum(2*((-1).^(1:n_cheb-1)).'.*cB(2:end,s));
    nn = (1:n_cheb-1).';
    Rs(2:n_cheb) = Rs(2:n_cheb) + 2*nn.*cB(2:end,s) + mu_s*cB(1:end-1,s);
    Rs(2:n_cheb-1) = Rs(2:n_cheb-1) - mu_s*cB(3:end,s);
    Rs(n_cheb+1) = Rs(n_cheb+1) + mu_s*cB(end,s);
    u = details.G{s}*Rs;
    PF(:,s) = u(1:n_cheb);
    dl = abs(mu_s)*details.rho1(s)*abs(u(end))/(1-details.rho(s));
    Y_3 = Y_3 + w_sp(s)*(sum(w(n_cheb+1:end).*abs(u(n_cheb+1:end))) + dl);
    dtot = dtot + w_sp(s)*dl;
end
A0_F = intval(details.A0)*PF(:);
omega = build_omega_vec_local(Mc);
Y_0 = sum(omega.*abs(A0_F)) + details.A0n*dtot;
Z_leN = intval(details.Z_leN); Z_A0 = intval(details.Z_A0); Z_inf = intval(details.Z_inf);
Y = sup(Y_0*(1 + Z_leN*(1+Z_inf)) + Y_3*(1+Z_inf)*(Z_A0 + 1));
Y_0 = sup(Y_0); Y_3 = sup(Y_3);
end

function w_sp = build_spatial_weights(Mc_sp)
    num_sp = prod(Mc_sp);
    sp_sz = Mc_sp;
    if isscalar(sp_sz), sp_sz = [sp_sz, 1]; end
    sp_subs = cell(1, max(length(Mc_sp), 1));
    [sp_subs{:}] = ind2sub(sp_sz, (1:num_sp)');
    w_sp = intval(ones(num_sp, 1));
    for d = 1:length(Mc_sp)
        for s = 1:num_sp
            if sp_subs{d}(s) > 1
                w_sp(s) = w_sp(s) * 2;
            end
        end
    end
end

% =========================================================================
% Helper: build omega vector (for Y_0)
% =========================================================================
function omega = build_omega_vec_local(Mc)
    n_cheb = Mc(end);
    Mc_sp = Mc(1:end-1);
    num_sp = prod(Mc_sp);
    total_dim = num_sp * n_cheb;
    dim_sp = length(Mc_sp);

    sp_sz = Mc_sp;
    if isscalar(sp_sz), sp_sz = [sp_sz, 1]; end
    sp_subs = cell(1, max(dim_sp, 1));
    [sp_subs{:}] = ind2sub(sp_sz, (1:num_sp)');

    w_sp = intval(ones(num_sp, 1));
    for d = 1:dim_sp
        for s = 1:num_sp
            if sp_subs{d}(s) > 1
                w_sp(s) = w_sp(s) * 2;
            end
        end
    end

    w_cheb = 2 * intval(ones(n_cheb, 1));
    w_cheb(1) = intval(1);

    omega = intval(zeros(total_dim, 1));
    for s = 1:num_sp
        for n = 1:n_cheb
            idx = (s-1)*n_cheb + n;
            omega(idx) = w_sp(s) * w_cheb(n);
        end
    end
end
