function [correction,pn,sigma_product] = getting_multistep_correction(Phi,sigma,remainder)
% Additive endpoint correction, independent of any inverse endpoint map.
% Phi holds SIGNED weighted finite enclosures. remainder(:,:,k) bounds the
% two-by-two operator block norms of U_k-E_k. Entries of correction{i}
% are two-by-two blocks for starts 1,...,i (packed horizontally).
K=numel(Phi);pn=intval(zeros(K));sigma_product=pn;
for hi=1:K
    product=intval(eye(size(Phi{hi},1)));tail=intval(1);
    for lo=hi:-1:1
        product=product*intval(Phi{lo});tail=tail*intval(sigma(lo));
        pn(hi,lo)=norm(product,1);sigma_product(hi,lo)=abs(tail);
    end
end
correction=cell(K,1);resolved=cell(K,1);
for hi=1:K
    row=intval(zeros(2,2*hi));
    for ell=1:hi
        if ell==hi,P=intval(eye(2));
        else,P=diag([pn(hi,ell+1),sigma_product(hi,ell+1)]);end
        R=P*intval(remainder(:,:,ell));
        row(:,2*ell-1:2*ell)=row(:,2*ell-1:2*ell)+R;
        if ell>1
            row(:,1:2*(ell-1))=row(:,1:2*(ell-1))+R*resolved{ell-1};
        end
    end
    correction{hi}=intval(max(0,sup(row)));row=correction{hi};
    for lo=1:hi
        row(:,2*lo-1:2*lo)=row(:,2*lo-1:2*lo) ...
            +diag([pn(hi,lo),sigma_product(hi,lo)]);
    end
    resolved{hi}=intval(max(0,sup(row)));
end
end
