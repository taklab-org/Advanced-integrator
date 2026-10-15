function psi = psi_k(ba_sup,k1,k2,m,nu) % ki starts from index 0
% size_ba = size(ba);
% n = size_ba(1); N1 = size_ba(2); N2 = size_ba(3);
[N1,N2] = size(ba_sup);
n1 = N1+k1; n2 = N2+k2;
%% Index ell
[ell1,ell2] = ndgrid(1:n1,1:n2); % indices F_{k+N} from one (because of picking off nonzero indices)
m1 = m(1); m2 = m(2);
ell1(1:m1,1:m2) = 0;
ell2(1:m1,1:m2) = 0;
ell1 = nonzeros(ell1)-1;
ell2 = nonzeros(ell2)-1; % F_{k+N}\F_m from zero up to n2-1
%
%% Index k+ell, k-ell
k1_minus_ell1 = abs(k1-ell1); k1_minus_ell1(k1_minus_ell1>N1) = N1; % up to N1-1. If index is N1, the elemtent must be zero
k1_plus_ell1 = abs(k1+ell1);  k1_plus_ell1(k1_plus_ell1 > N1) = N1;
k2_minus_ell2 = abs(k2-ell2); k2_minus_ell2(k2_minus_ell2>N2) = N2;
k2_plus_ell2 = abs(k2+ell2);  k2_plus_ell2(k2_plus_ell2 > N2) = N2;
%
%% weight
% alp = 4*ones(N1+1,N2+1); alp(1,:) = 2; alp(:,1) = 2; alp(1,1) = 1;
% v1 = nu.^(0:N1); v2 = nu.^(0:N2);
% w = (v1'*v2) .* alp;
% w_vec = w(sub2ind([N1+1,N2+1],ell1+1,ell2+1));
alp = 4*ones(n1,n2); alp(1,:) = 2; alp(:,1) = 2; alp(1,1) = 1;
v1 = nu.^(0:n1-1); v2 = nu.^(0:n2-1);
% w = (v1'*v2) .* alp;
w = alp .* (v1'*v2);
w_vec = w(sub2ind([n1,n2],ell1+1,ell2+1));
%
% %========= computing with chebyshev
% tic
% ba2_ext = intval(zeros(n+1,N1+1,N2+1));
% ba2_ext(1:n,1:N1,1:N2) = ba2;
% %
% nind = kron((1:n)',ones(length(ell1),1));
% N1ind_minus = kron(ones(n,1),k1_minus_ell1+1);
% N2ind_minus = kron(ones(n,1),k2_minus_ell2+1);
% ba2_k_minus_ell = reshape(ba2_ext(sub2ind(size(ba2_ext),nind,N1ind_minus,N2ind_minus)),length(ell1),n);
% %
% N1ind_plus = kron(ones(n,1),k1_plus_ell1+1);
% N2ind_plus = kron(ones(n,1),k2_plus_ell2+1);
% ba2_k_plus_ell = reshape(ba2_ext(sub2ind(size(ba2_ext),nind,N1ind_plus,N2ind_plus)),length(ell1),n);
% %
% psi = max(mag(chebmag((ba2_k_minus_ell + ba2_k_plus_ell)')./w_vec))
% toc
% %===================
%
% tic
% ba_sup = reshape(chebmag(reshape(ba,n,N1*N2)),N1,N2);
ba_sup_ext = intval(zeros(N1+1,N2+1));
ba_sup_ext(1:N1,1:N2) = ba_sup;
% ba_sup_k_minus_ell = ba_sup_ext(sub2ind([N1+1,N2+1],k1_minus_ell1+1,k2_minus_ell2+1));
% ba_sup_k_plus_ell  = ba_sup_ext(sub2ind([N1+1,N2+1],k1_plus_ell1+1,k2_plus_ell2+1));
% psi = max(mag(max(ba_sup_k_minus_ell, ba_sup_k_plus_ell)./w_vec));
% toc
tmp = max(max(ba_sup_ext(sub2ind([N1+1,N2+1],k1_minus_ell1+1,k2_minus_ell2+1)),ba_sup_ext(sub2ind([N1+1,N2+1],k1_plus_ell1+1,k2_minus_ell2+1))),...
  max(ba_sup_ext(sub2ind([N1+1,N2+1],k1_minus_ell1+1,k2_plus_ell2+1)) + ba_sup_ext(sub2ind([N1+1,N2+1],k1_plus_ell1+1,k2_plus_ell2+1))));
psi = max(mag(tmp./w_vec));
