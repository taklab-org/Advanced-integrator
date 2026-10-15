function cc = pdecheb2(f, I, tol, Nmax)% return Two-sided Chebyshev
% f: ode function (Fourier coefficients) by ode45, ode78, ode89 etc.
% I: domain of the function f
% tol: tolerance of Chebyshev coefficients
% Nmax: maximum # of Chebyshev coefficients
arguments
    f
    I = [-1,1];
    tol = 1e-15;
    Nmax = 100;
end
%
devalFun = @(x) deval(f, x).';
%
i = 3;
%
while true
    test_size = 2^i+1;
    schbc = chebcoeffs(devalFun,test_size,I);
    if norm(schbc(end-4:end,:),inf) < tol || (test_size > Nmax)
      [nonzero_row_index,~]=find(abs(schbc)>tol);
      chopsize = max(nonzero_row_index);
      break
    end
    % [~,nonzero_row_index]=find(abs(schbc)>tol);
    % sample_column = schbc(:,min(nonzero_row_index));
    % chopsize = standardChop(sample_column,tol/abs(sample_column(1)));
    % if test_size > chopsize || (test_size > Nmax)
    %     break
    % end
    i = i + 1;
end
cc = schbc(1:min(Nmax,chopsize+1),:);
end
