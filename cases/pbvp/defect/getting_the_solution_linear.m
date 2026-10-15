function [ba1,ba2,ba3,ba4,delta,a_end] = getting_the_solution_linear(ba,N,tspan,a0,f_lin1,f_lin2,f_lin3,f_lin4,params,L,nu)
%%
lambda  = params(1);
N1 = N(1); N2 = N(2);
h = tspan(2)-tspan(1);
%
% a0 = reshape(a0_init,N1*N2,1);
%
% addpath("defect/Chebyshev/")
% opts = odeset('abstol',1e-16,'reltol',2.22045e-14);
% f = ode15s(@(t,y) F_SH_1D(y,params,N,mid(L)),tspan,a0,opts);
% f = ode113(@(t,y) F_SH_2D(y,lambda,[N1,N2]),tspan,a0,opts);
a1 = pdecheb2(f_lin1,tspan); % Two-sided Chebyshev
a2 = pdecheb2(f_lin2,tspan); % Two-sided Chebyshev
a3 = pdecheb2(f_lin3,tspan); % Two-sided Chebyshev
a4 = pdecheb2(f_lin4,tspan); % Two-sided Chebyshev
n1 = size(a1,1); n2 = size(a2,1);  n3 = size(a3,1);  n4 = size(a4,1);
disp(['(n1, n2, n3, n4) = ',num2str([n1,n2,n3,n4])])
%
Chebyshev = size(ba,1); Fourier = prod(N);
% if size_ba_cheb ~= n1
[ba1,a1,n1] = padding_ba_a(reshape(ba,[Chebyshev,Fourier]),a1);
a1_full = (zeros(n1,N1*N2)); a1_full(1,1) = lambda; a1_full(:,2:end) = mid(a1);
[ba2,a2,n2] = padding_ba_a(reshape(ba,[Chebyshev,Fourier]),a2);
a2_full = (zeros(n2,N1*N2)); a2_full(1,1) = lambda; a2_full(:,2:end) = mid(a2);
[ba3,a3,n3] = padding_ba_a(reshape(ba,[Chebyshev,Fourier]),a3);
a3_full = (zeros(n3,N1*N2)); a3_full(1,1) = lambda; a3_full(:,2:end) = mid(a3);
[ba4,a4,n4] = padding_ba_a(reshape(ba,[Chebyshev,Fourier]),a4);
a4_full = (zeros(n4,N1*N2)); a4_full(1,1) = lambda; a4_full(:,2:end) = mid(a4);
%
%
%%%%%%% Refine a via Newton: OUTPUT is Two-sided Chebyshev
[a1, a1_end] = refine_newton_linear(a1_full.',mid(ba1).',a0(:,1),n1,N,L,params,h); % size(a) = [N, n], size(a0) = [N, 2]
a1 = a1.';
[a2, a2_end] = refine_newton_linear(a2_full.',mid(ba2).',a0(:,2),n2,N,L,params,h); % size(a) = [N, n], size(a0) = [N, 2]
a2 = a2.';
[a3, a3_end] = refine_newton_linear(a3_full.',mid(ba3).',a0(:,3),n3,N,L,params,h); % size(a) = [N, n], size(a0) = [N, 2]
a3 = a3.';
[a4, a4_end] = refine_newton_linear(a4_full.',mid(ba4).',a0(:,4),n4,N,L,params,h); % size(a) = [N, n], size(a0) = [N, 2]
a4 = a4.';
%TODO: newton does not converge once....???

a_end = [a1_end, a2_end, a3_end, a4_end];
%
a1 = intval(a1);
a2 = intval(a2);
a3 = intval(a3);
a4 = intval(a4);
% plotcoeffs(u_cheb),return
% mesh(abs(reshape(u_cheb(end),N1,N2)))
% set(gca, 'ZScale', 'log')
% return
%%%%%%% Compute residual
% ta = [a(1,:);a(2:end,:)/2];
%
% time derivative
rescaleFactork = h/2;
du1 = chebdiff(a1)/rescaleFactork;
du1 = [du1(1,:);du1(2:end,:)/2];
du2 = chebdiff(a2)/rescaleFactork;
du2 = [du2(1,:);du2(2:end,:)/2];
du3 = chebdiff(a3)/rescaleFactork;
du3 = [du3(1,:);du3(2:end,:)/2];
du4 = chebdiff(a4)/rescaleFactork;
du4 = [du4(1,:);du4(2:end,:)/2];

%
if exist('intval','file') && isintval(a1(1))
  du1_ext = intval(zeros(3*n1-2,3*N1-2,3*N2-2));
  du2_ext = intval(zeros(3*n2-2,3*N1-2,3*N2-2));
  du3_ext = intval(zeros(3*n3-2,3*N1-2,3*N2-2));
  du4_ext = intval(zeros(3*n4-2,3*N1-2,3*N2-2));
else
  du1_ext = zeros(3*n1-2,3*N1-2,3*N2-2);
  du2_ext = zeros(3*n2-2,3*N1-2,3*N2-2);
  du3_ext = zeros(3*n3-2,3*N1-2,3*N2-2);
  du4_ext = zeros(3*n4-2,3*N1-2,3*N2-2);
end
du1_ext(1:n1-1,1:N1,1:N2) = reshape(du1,n1-1,N1,N2);
du2_ext(1:n2-1,1:N1,1:N2) = reshape(du2,n2-1,N1,N2);
du3_ext(1:n3-1,1:N1,1:N2) = reshape(du3,n3-1,N1,N2);
du4_ext(1:n4-1,1:N1,1:N2) = reshape(du4,n4-1,N1,N2);
%
%
ta1 = [a1(1,:);a1(2:end,:)/2]; % Two-sided --> One-sided Chebyshev
Fa1_ext = F_OK_2D_linear(ta1,ba1,params,n1,N,L);
ta2 = [a2(1,:);a2(2:end,:)/2]; % Two-sided --> One-sided Chebyshev
Fa2_ext = F_OK_2D_linear(ta2,ba2,params,n2,N,L);
ta3 = [a3(1,:);a3(2:end,:)/2]; % Two-sided --> One-sided Chebyshev
Fa3_ext = F_OK_2D_linear(ta3,ba3,params,n3,N,L);
ta4 = [a4(1,:);a4(2:end,:)/2]; % Two-sided --> One-sided Chebyshev
Fa4_ext = F_OK_2D_linear(ta4,ba4,params,n4,N,L);
%
%
defect1 = du1_ext - Fa1_ext;
defect2 = du2_ext - Fa2_ext;
defect3 = du3_ext - Fa3_ext;
defect4 = du4_ext - Fa4_ext;
delta = [wnorm(defect1,nu), wnorm(defect2,nu), wnorm(defect3,nu), wnorm(defect4,nu)];
%
ba1 = reshape(ta1,n1,N1,N2);
ba2 = reshape(ta2,n2,N1,N2);
ba3 = reshape(ta3,n3,N1,N2);
ba4 = reshape(ta4,n4,N1,N2);
%
end
