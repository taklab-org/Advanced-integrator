function [F] = F_lin_fin_dim(mu_k,Nop_coeff,a_bar,c,h,e_j,L,lambda)

if ~exist('lambda','var')
    lambda = [];
end

Mc = size(c);
Ma = size(a_bar);
c_temp = c;
mu_k_temp = mu_k;
c = zeros(Ma);
mu_k = zeros(Ma);
for i = 1:length(Ma)
    otherdims_a{i} = 1:Ma(i);
    otherdims_c{i} = 1:Mc(i);
    otherdims_zeros{i} = 1;
end
c(otherdims_c{:}) = c_temp;
mu_k(otherdims_c{:}) = mu_k_temp;
N = size(a_bar)-1;

a_bar_pad = zeros(Ma+1);
a_bar_pad(otherdims_a{:}) = a_bar;

c_pad = zeros(Ma+1);
c_pad(otherdims_a{:}) = c;

mu_k_pad= zeros(Ma+1);
mu_k_pad(otherdims_a{:}) = mu_k;

% This may be the same code as above
% asize = size(a_bar);
% pad_size = 1;
% a_bar_pad = zeros(asize+pad_size);
% a_bar_pad(1:end-pad_size,1:end-pad_size,1:end-pad_size,1:end-pad_size) = a_bar;
% c_pad = zeros(asize+pad_size);
% c_pad(1:end-pad_size,1:end-pad_size,1:end-pad_size,1:end-pad_size) = c;
% mu_k_pad = zeros(asize+pad_size);
% mu_k_pad(1:end-pad_size,1:end-pad_size,1:end-pad_size,1:end-pad_size) = mu_k;


nop = N_op_lin(Nop_coeff,a_bar_pad,c_pad);
M = size(nop)-1;
dim = length(M)-1;
s = arrayfun(@(k) 0:k, M(1:end-1), 'UniformOutput', false);
k = cell(1, dim);
[k{:}] = ndgrid(s{:});
if (size(Nop_coeff,1) - 1 ) >0
    K_sq = zeros(size(k{1}));
    for ell = 1:dim
        K_sq = K_sq + (k{ell}*L(ell)).^(size(Nop_coeff,1) - 1 );
    end
else
    K_sq = ones(size(k{1}));
end
N_rep = ones(size(M));
N_rep(end) = M(end)+1;
K_sq = repmat(K_sq,N_rep);
nop = K_sq.*nop;

%
% indstring=repmat('1:end-2,',[1 ndims(nop)]);
% indstring(end)=[];
% eval(['nop=nop(' indstring ');']);




psi = h/2*(mu_k_pad.*c_pad + nop);
F = lambda_op(c_pad) + T_op(psi);
F = F(otherdims_c{:});

% indstring=repmat('1:end-1,',[1 ndims(nop)]);
% indstring(end)=[];
% eval(['F=F(' indstring ');']);

% otherdims = repmat({':'},1,length(N)-1);
F(otherdims_c{1:(end-1)}, 1) =  - e_j;
for j = 0:N(end)
    if j == 0
        F(otherdims_c{1:end-1}, 1) = F(otherdims_c{1:end-1}, 1) + c(otherdims_c{1:end-1}, j+1);
    else
        F(otherdims_c{1:end-1}, 1) = F(otherdims_c{1:end-1}, 1) + 2*((-1)^j)*c(otherdims_c{1:end-1}, j+1);
    end
    if ~isempty(lambda) && j > 0
        F(otherdims_zeros{1:end-1},j+1) = 2*j*c(otherdims_zeros{1:end-1},j+1) ;
    end

end


end
