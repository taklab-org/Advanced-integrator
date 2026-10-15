function [F] = F_lin_fin_dim_int(mu_k,Nop_coeff,a_bar,c,h,e_j,L,lambda)
if ~exist('lambda','var')
    lambda = [];
end
Mc = size(c);
Ma = size(a_bar);
c_temp = c;
mu_k_temp = mu_k;
c = intval(zeros(Ma));
mu_k = intval(zeros(Ma));
for i = 1:length(Ma)
    otherdims_a{i} = 1:Ma(i);
    otherdims_c{i} = 1:Mc(i);
    otherdims_zeros{i} = 1;
end
c(otherdims_c{:}) = c_temp;
mu_k(otherdims_c{:}) = mu_k_temp;
N = size(a_bar)-1;

a_bar_pad = intval(zeros(Ma+1));
a_bar_pad(otherdims_a{:}) = a_bar;

c_pad = intval(zeros(Ma+1));
c_pad(otherdims_a{:}) = c;

mu_k_pad= intval(zeros(Ma+1));
mu_k_pad(otherdims_a{:}) = mu_k;

% a_bar_pad = a_bar;
% a_bar_pad(otherdims_ext{:},N(end)+2) = 0 ;
% c_pad = c;
% c_pad(otherdims_ext{:},N(end)+2) = 0 ;
% mu_k_pad = mu_k;
% mu_k_pad(otherdims_ext{:},N(end)+2) = 0 ;

nop = N_op_lin_int(Nop_coeff,a_bar_pad,c_pad);
M = size(nop)-1;
dim = length(M)-1;
s = arrayfun(@(k) 0:k, M(1:end-1), 'UniformOutput', false);
k = cell(1, dim);
[k{:}] = ndgrid(s{:});
if (size(Nop_coeff,1) - 1 ) >0
    K_sq = intval(zeros(size(k{1})));
for ell = 1:dim
    K_sq = K_sq + (k{ell}*L(ell)).^(size(Nop_coeff,1) - 1 );
end
else
    K_sq = intval(ones(size(k{1})));
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
F = lambda_op_int(c_pad) + T_op_int(psi);
F = F(otherdims_c{:});

% indstring=repmat('1:end-1,',[1 ndims(nop)]);
% indstring(end)=[];
% eval(['F=F(' indstring ');']);
otherdims = repmat({':'},1,length(N)-1);
% [Opt7] Boundary condition: separate j=0 case and lambda loop
F(otherdims_c{1:end-1}, 1) = F(otherdims_c{1:end-1}, 1) - e_j + c(otherdims_c{1:end-1}, 1);
for j = 1:N(end)
    F(otherdims_c{1:end-1}, 1) = F(otherdims_c{1:end-1}, 1) + 2*((-1)^j)*c(otherdims_c{1:end-1}, j+1);
end
if ~isempty(lambda)
    for j = 1:N(end)
        F(otherdims_zeros{1:end-1},j+1) = 2*j*c(otherdims_zeros{1:end-1},j+1);
    end
end
end
