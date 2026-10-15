function DF = DF_map(a,n,N,L,params,h)
% Input: size(a)  = [n*N1*N2, 1] : Vector form (One-sided Chebyshev*cosFourier*cosFourier, 1)
% Output: size(DF) = [n*N1*N2, n*N1*N2] : Jacobian Matrix
%
%
% lambda  = params(1);
sigma   = params(2);
epsilon = params(3);
% Wave numbers: N
N1 = N(1); N2 = N(2);
N1N2 = prod(N);
L1 = L(1); L2 = L(2);
%
% Convert Vector form -> Matrix form: (cosFourier*cosFourier, One-sided Chebyshev)
a = reshape(a,N1N2,n);
%
%
%% Linear part
%%%%%%%%%%%%%%%%%
% Time derivative
%%%%%%%%%%%%%%%%%
Dt_mat = chebdiff_oneside(eye(n),h);
%
% Derivative for initial value
ell = 0:n-1; v = (-1).^ell; v(2:end) = 2*v(2:end);
%
% Linear operator (time)
Dt = kron([Dt_mat;v],eye(N1N2));
% Dt = kron(eye(N1N2N3),[Dt_mat;v]);
%
%%%%%%%%%%%%%%%%%
% Sace derivative
%%%%%%%%%%%%%%%%%
% [n1, n2, n3] = ndgrid((-N1+1):(N1-1),(-N2+1):(N2-1),(-N3+1):(N3-1));
% k = (0:N-1)';
[kx,ky] = ndgrid(0:N1-1,0:N2-1);
%
bkL2 = (kx*L1).^2 + (ky*L2).^2;
%
% Laplacian operator
D2 = [repmat(bkL2(:),[n-1,1]); zeros(N1N2,1)]; % vector form [n*N1N2, 1]
% Linear operator
% mu_k = -beta1*(k*L).^4 - beta2*(k*L).^2 + 1;
% mu_k = lambda - (1-(kx*L1).^2-(ky*L2).^2).^2 ;
mu_k = bkL2 .* (-epsilon^2*bkL2 + 1) - sigma;
%
%
L = diag([repmat(mu_k(:),[n-1,1]); zeros(N1N2,1)]); % diagonal operator
% L  = permute(repmat(mu_k,[1,1,1,n]),[4,1,2,3]);  L = diag(L(:));
%
%% Convolutions
% Convert Matrix form: (cosFourier*cosFourier, One-sided Chebyshev) -> Tensor form: (cosFourier, cosFourier, One-sided Chebyshev)
a_tensor = reshape(a,[N1,N2,n]);
a2 = powerconvcos(a_tensor,2); % Tensor form: (cosFourier, cosFourier, One-sided Chebyshev)
DFa = -3*a2;
% DFa = 3*[a2(2:end,:);a2]; % 3 * bar{a} & Convert to (Fourier, One-sided Chebyshev)
FoChTa = [0 3 0]; Nout = [N-1,n-1];
[~,DRHS] = multiplication(DFa, FoChTa, Nout);
%
%%%%%%%%%%%%%%%%%
% DF
%%%%%%%%%%%%%%%%%
DF = Dt - (L + D2.*DRHS);
end
