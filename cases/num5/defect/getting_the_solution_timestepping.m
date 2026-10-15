function [ba, delta, a_end] = getting_the_solution_timestepping(N,tspan,a0_init,params,L,nu,rigorous)
% params = [lambda,sigma,epsilon];
lambda  = params(1);
% sigma   = params(2);
% epsilon = params(3);
%%
N1 = N(1); N2 = N(2);
h = tspan(2)-tspan(1);
%
a0 = reshape(a0_init, N1*N2,1);
a0 = a0(2:end); % remove the zero mode
%
%
% addpath("defect/Chebyshev/")
opts = odeset('abstol',1e-16,'reltol',2.22045e-14);
% opts = odeset('abstol',1e-11,'reltol',2.22045e-14);
% f = ode15s(@(t,y) F_DC_2D(y,params,N,L),tspan,a0,opts); % f is supposed to except zero mode
% f = ode78(@(t,y) F_DC_2D(y,params,N,L),tspan,a0,opts); % f is supposed to except zero mode
f = ode113(@(t,y) F_OK_2D(y,params,N,L),tspan,a0,opts);
a = pdecheb2(f,tspan); % Two-sided Chebyshev
n = size(a,1); disp(['n = ',num2str(n)])
a_end = reshape([lambda,f.y(:,end).'],N1,N2); % a_end includes the zero mode
% 2.9422e-09(113), 2.9405e-09(89), 3.2228e-09(78), 2.9404e-09(45)
%
if rigorous
  a_full = intval(zeros(n,N1*N2));
  a_full(1,1) = lambda;
  a_full(:,2:end) = a; % a_full includes the zero mode
else
  a_full = zeros(n,N1*N2);
  a_full(1,1) = lambda;
  a_full(:,2:end) = a; % a_full includes the zero mode
end
%
%%%%%%% Compute residual
% time derivative
rescaleFactork = h/2;
du = chebdiff(a_full)/rescaleFactork;
du = [du(1,:);du(2:end,:)/2];
du = reshape(du,n-1,N1,N2); % du includes the zero mode
%
if exist('intval','file') && isintval(a_full(1))
  du_ext = intval(zeros(3*n-2,3*N1-2,3*N2-2));
else
  du_ext = zeros(3*n-2,3*N1-2,3*N2-2);
end
du_ext(1:n-1,1:N1,1:N2) = du;
%
ta = [a_full(1,:);a_full(2:end,:)/2]; % Two-sided --> One-sided Chebyshev
Fa_ext = F_OK_2D_ext(ta,params,n,N,L);
%
defect = du_ext - Fa_ext;
defect(2:end,:,:) = 2*defect(2:end,:,:); % back to Two-sided Chebyshev
defect(:,1,1) = 0; % defect of zero mode must be 0 because of the convervation
delta = wnorm(permute(sum(abs(defect),1),[2,3,1]),nu);
% delta = 0;
ba = reshape(ta,n,N1,N2); % return One-sided Chebyshev including the zero mode
if rigorous
    ba = mid(ba);
end
end
