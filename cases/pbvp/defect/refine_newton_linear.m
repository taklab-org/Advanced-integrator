function [x, a_end] = refine_newton_linear(x,ba,x0,n,N,L,params,h)
% Input: size(x)  = [N, n] : Matrix form (cosFourier, Two-sided Chebyshev)
% Output: size(x) = [N, n] : Matrix form (cosFourier, Two-sided Chebyshev)
%
% tol = 5e-15;
tol = eps*numel(x);
%
% Two-sided --> One-sided Chebyshev
x(:,2:end) = x(:,2:end)/2;
%
% Convert Matrix form -> Vector form
[Fourier, Chebyshev] = size(x);
%
x = x(:);
ba_x = ba(:); % One-sided Chebyshev
% % Construct Full Fourier
% x_fullFourier = [flipud(x(2:end,:));x];
% x = x_fullFourier(:);
%
% The zero map F (cosFourier)
F = F_map_linear(x,ba,x0,n,N,L,params,h); % One-sided Chebyshev -> One-sided Chebyshev
%
% display(['Before step #1, ||F||_1 = ',num2str(norm(F,1))])
itercount = 0;
%
while (itercount<=10) && (norm(F,1) > tol)
    % DF_diff = finite_diff_DF(x,x0,n,N,L,params,h);
    DF = DF_map(ba_x,n,N,L,params,h);
    % display(['Size of DF: ',num2str(size(DF))])
    x  = x - DF\F;
    F  = F_map_linear(x,ba,x0,n,N,L,params,h);
    display(['Newton step # ',num2str(itercount+1),', ||F||_1 = ',num2str(norm(F,1)),', ||DF||_1 = ',num2str(norm(DF,1))])
    itercount = itercount+1;
end
%
% Convert Vector form -> Matrix form
x = reshape(x,Fourier,Chebyshev);
% x = reshape(permute(reshape(x, N1, n, 3), [1,3,2]), Fourier,Chebyshev);
x(:,2:end) = 2*x(:,2:end); % back to Two-sided Chebyshev
%
a_end = sum(x,2);
