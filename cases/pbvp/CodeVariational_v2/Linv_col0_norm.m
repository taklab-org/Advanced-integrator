function nrm = Linv_col0_norm(N, mu)
% LINV_COL0_NORM  Compute ||W_0 (L_k^{(N)})^{-1}_{col(0)}||_1
%
%   This is the weighted 1-norm of the first column of L_k^{(N)} inverse,
%   needed for the Y_2 bound (eq. 34 in the paper).
%   W_0 is the diagonal operator with W_0(v_0) = v_0, W_0(v_n) = 2*v_n for n>=1.
%
%   For the finite-dimensional case, we compute this directly from the
%   matrix inverse.

L = complete_op(N, mu);
Linv = inv(L);
col0 = Linv(:, 1);   % first column

% Apply W_0 weights: 1 for n=0, 2 for n>=1
w = 2 * ones(N+1, 1);
if isa(mu, 'intval')
    w = intval(w);
end
w(1) = 1;

nrm = sum(abs(col0) .* w);

end
