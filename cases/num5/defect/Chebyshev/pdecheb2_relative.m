function cc = pdecheb2_relative(f, I, rel_tol, abs_tol, Nmax)
% pdecheb2 with per-mode relative tolerance criterion.
% Stop when, for every column k of schbc, max|tail_k| / max|column_k| < rel_tol.
% Also enforce absolute tail < abs_tol as a safety net.
%
% f: ode result (deval-able)
% I: domain
% rel_tol: relative tolerance per mode (default 1e-10)
% abs_tol: absolute floor (default 1e-15) -- if max amplitude is below this,
%          treat that mode as "noise" and ignore in stop criterion
% Nmax: maximum # of Chebyshev coefficients
arguments
    f
    I = [-1,1];
    rel_tol = 1e-10;
    abs_tol = 1e-15;
    Nmax = 100;
end

devalFun = @(x) deval(f, x).';

i = 3;
while true
    test_size = 2^i + 1;
    schbc = chebcoeffs(devalFun, test_size, I);
    amp_k = max(abs(schbc), [], 1);          % 1xN row vector of per-mode peak
    tail_k = max(abs(schbc(end-4:end, :)), [], 1);
    % Only check modes whose amplitude exceeds the noise floor
    active = amp_k > abs_tol;
    if all(~active)
      relmax = 0;
    else
      relmax = max(tail_k(active) ./ amp_k(active));
    end
    if relmax < rel_tol || (test_size > Nmax)
      [nonzero_row_index, ~] = find(abs(schbc) > abs_tol);
      if isempty(nonzero_row_index)
        chopsize = 1;
      else
        chopsize = max(nonzero_row_index);
      end
      break
    end
    i = i + 1;
end

cc = schbc(1:min(Nmax, chopsize+1), :);
end
