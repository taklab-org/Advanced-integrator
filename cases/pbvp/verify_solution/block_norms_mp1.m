function [u11, u12, u21, u22] = block_norms_mp1(A, m)
% Extract the block norms from an (m+1)x(m+1) interval matrix.
% Returns sup (double) for fast storage in Step 2.
% Converted back to intval in bulk before Step 3.
u11 = sup(max(sum(abs(A(1:m,1:m)), 1)));
u12 = sup(sum(abs(A(1:m,m+1))));
u21 = sup(max(abs(A(m+1,1:m))));
u22 = sup(abs(A(m+1,m+1)));
end
