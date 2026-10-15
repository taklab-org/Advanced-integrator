function plot_cheb_fourier(a,t_int,x_int,index)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Input is the two-sided Chebyshev     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a: (Chebyshev, Fourier) array
% t_int: time interval
% x_int: space interval
% index 1: plot profile
%       2: plot Fourier modes
%       3: plot maximum norm
%       4: plot Chebyshev coeffs
%
% SPACE evaluation
[M,n] = size(a); % n=2N-1
N = (n+1)/2;
%
% TIME evaluation
ta = t_int(1); tb = t_int(2);
%
%
if index==1
    % Plot profile:
    xa = x_int(1); xb = x_int(2);
    n_pad = 200;
    a_pad = [zeros(M,n_pad), a, zeros(M,n_pad)];
    N_pad = N + n_pad;
    h_pad = (xb-xa)/(2*N_pad-1);
    x_pad = xa + h_pad*(0:2*N_pad-2);
    ax_pad = (2*N_pad-1)*ifft(ifftshift(a_pad,2).').';
    t = linspace(ta,tb,1000)';
    axt = eval_cheb(ax_pad,t);
    %
    subplot(1,2,1);
    mesh(x_pad,t,real(axt),'EdgeColor','none','FaceColor','flat')
    %
    xlabel('$x$','interpreter','latex'), ylabel('$t$','interpreter', 'latex'), zlabel('$\mathrm{Re}(u)$','interpreter', 'latex')
    title('Real part')
    %
    subplot(1,2,2);
    mesh(x_pad,t,imag(axt),'EdgeColor','none','FaceColor','flat')
    %
    xlabel('$x$','interpreter','latex'), ylabel('$t$','interpreter', 'latex'), zlabel('$\mathrm{Im}(u)$','interpreter', 'latex')
    title('Imaginary part')
    %
elseif index==2
    % Plot fourier modes:
    t = linspace(ta,tb,100)';
    k = ((-N+1):(N-1))';
    at = eval_cheb(a,t);
    subplot(1,2,1);
    mesh(k,t,abs(real(at)))
    %
    xlabel('$k$','interpreter','latex'), ylabel('$t$','interpreter', 'latex'), zlabel('$\mathrm{Re}(u)$','interpreter', 'latex')
    title('Real part')
    set(gca, 'ZScale', 'log')
    %
    subplot(1,2,2);
    mesh(k,t,abs(imag(at)))
    %
    xlabel('$k$','interpreter','latex'), ylabel('$t$','interpreter', 'latex'), zlabel('$\mathrm{Im}(u)$','interpreter', 'latex')
    title('Imaginary part')
    set(gca, 'ZScale', 'log')
    %
elseif index==3
    % maximum norm plot
    t = linspace(ta,tb,100)';
    LW = 'linewidth'; lw = 1.6;
    y=ifft(ifftshift(eval_cheb(a,t),2).');
    plot(t,abs(max(y)),LW,lw);
    xlabel('$t$','interpreter','latex'), ylabel('$\|u(t)\|_{\infty}$','interpreter', 'latex')
    title('The maximum norm of the solution.')
    %
elseif index==4
    % plot Chebyshev coeffs
    plot_chebcoeffs(a)
    title('Chebyshev coefficients.')
    %
end
