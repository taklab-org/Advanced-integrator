function [ba, defect, a_end] = getting_the_solution_refine_newton(N,tspan,a0,params,L)
% params = [lambda,sigma,epsilon];
lambda  = params(1);
% sigma   = params(2);
% epsilon = params(3);
%%
N1 = N(1); N2 = N(2);
h = tspan(2)-tspan(1);
%
% disp('Compute numerical data ....')
tol = 1e-16;
opts = odeset('abstol',tol,'reltol',2.22045e-14);
% opts = odeset('abstol',tol,'reltol',tol);
% tic
Y = ode113(@(t,y) F_OK_2D(y,params,N,mid(L)),tspan,a0(2:end),opts);
% Y = ode15s(@(t,y) F_OK_2D(y,params,N,mid(L)),tspan,a0(2:end),opts);
% toc
%
a = pdecheb2(Y,tspan); % size(a) = [n, N1, N2]: Two-sided Chebyshev
n = size(a,1); disp(['# of Chebyshev: ',num2str(n)])
%
a_full = (zeros(n,N1*N2));
a_full(1,1) = lambda;
a_full(:,2:end) = a;
%
%
%%%%%%% Refine a via Newton: OUTPUT is Two-sided Chebyshev
[a_full, a_end] = refine_newton(a_full.',a0,n,N,L,params,h); % size(a) = [N, n], size(a0) = [N, 1]
a_full = a_full.';
%
a = intval(a_full);
% plotcoeffs(u_cheb),return
% mesh(abs(reshape(u_cheb(end),N1,N2)))
% set(gca, 'ZScale', 'log')
% return
%%%%%%% Compute residual
% ta = [a(1,:);a(2:end,:)/2];
%
% time derivative
rescaleFactork = h/2;
du = chebdiff(a)/rescaleFactork;
du = [du(1,:);du(2:end,:)/2];
du = reshape(du,n-1,N1,N2);
%
if exist('intval','file') && isintval(a(1))
  du_ext = intval(zeros(3*n-2,3*N1-2,3*N2-2));
else
  du_ext = zeros(3*n-2,3*N1-2,3*N2-2);
end
du_ext(1:n-1,1:N1,1:N2) = du;
%
%
ta = [a(1,:);a(2:end,:)/2]; % Two-sided --> One-sided Chebyshev
Fa_ext = F_OK_2D_ext(ta,params,n,N,L);
%
defect = du_ext - Fa_ext;
% delta = wnorm(defect,nu);
% delta = wnorm(intval(sup(abs(defect))),nu);
% defect(2:end,:) = 2*defect(2:end,:); % back to Two-sided Chebyshev
% delta = wnorm(sum(abs(defect),1),nu);
% delta3 = wnorm(reshape(chebmag(reshape(defect,3*n-2,(3*N1-2)*(3*N2-2)),tspan),(3*N1-2),(3*N2-2)),1)
%
% ba = ta;
ba = reshape(ta,n,N1,N2); % return One-sided Chebyshev
% rmpath("defect/Chebyshev/");
end
