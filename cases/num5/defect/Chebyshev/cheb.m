function cc = cheb(f, I, tol, Nmax)% return Two-sided Chebyshev
% f: given function handle
% I: domain of the function f
% tol: tolerance of Chebyshev coefficients
% Nmax: maximum # of Chebyshev coefficients
arguments
    f
    I = [-1,1];
    tol = eps;
    Nmax = 10000;
end
%
a = I(1); b = I(2);
% m = 0.5*(a+b); r = 0.5*(b-a); x = rand(5,1);
% x1 = m + x*r; x2 = m - x*r;
% fx1 = f(x1); fx2 = f(x2);
% For even/odd functions
f_sample = f(linspace(a,b,18)');
isEven = max(abs(f_sample-flipud(f_sample)),[],1) == 0;
isOdd = max(abs(f_sample+flipud(f_sample)),[],1) == 0;

%
% [~,num_of_fun] = size(f(a));
% odd_even = zeros(num_of_fun,1);
%
% for ell = 1:num_of_fun
%     if all(abs(fx1(:,ell)-fx2(:,ell))./abs(fx1(:,ell)) < 5*eps) % tolerance for odd/even function
%         odd_even(ell) = 1; % even function: 1
%     elseif all(abs(fx1(:,ell)+fx2(:,ell))./abs(fx1(:,ell)) < 5*eps)
%         odd_even(ell) = -1; %  odd function: -1
%     else
%         odd_even(ell) = 0; % otherwise: 0
%     end
% end
i = 3;
% schbc = 0; % sampling chebyshev coefficients
%
while true
    test_size = 2^i+1;
    schbc = chebcoeffs(f,test_size,I);
    chopsize = standardChop(schbc,tol);
    if test_size > chopsize || (2^i+1 > Nmax)
        break
    end
    i = i + 1;
end
% [M, ~] = find(abs(schbc) > tol);
% M = max(M);
% cc = chebcoeffs(f,M,I);
cc = schbc(1:chopsize,:);
% for ell = 1:num_of_fun
%     if odd_even(ell) == 1 % even function
%         cc(2:2:end,ell) = 0;
%     elseif odd_even(ell) == -1 % odd function
%         cc(1:2:end,ell) = 0;
%     end
% end
%
% adjust coefficients for symmetry
cc(2:2:end,isEven) = 0;
cc(1:2:end,isOdd) = 0;
end
