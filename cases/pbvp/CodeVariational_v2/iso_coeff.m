function [k] = iso_coeff(M,N)
dim =length(N);
k = zeros(1,dim);
for i = 1:dim
    k(i) = mod(M-1,N(i)+1);
    M = (M-k(i)-1)/(N(i)+1)+1;
end
end
