function [M] = sym_dif(k)
dim = length(k);
M = [k(1);-k(1)];
for i = 2:dim
    Mp = [M,k(i)*ones(size(M,1),1)];
    Mm = [M,-k(i)*ones(size(M,1),1)];
    M = [Mp;Mm];
end
M = unique(M,'rows');
end
