function cc = pdecheb2_relative_v2(f, I, rel_tol, Nmax)
% pdecheb2 with noise floor scaled by global max amplitude.
% Active mode: amp_k > rel_tol * max_global_amp.
% Stop when, for every active mode k, max|tail_k| / max|column_k| < rel_tol.
%
% f: ode result (deval-able)
% I: domain
% rel_tol: relative tolerance per mode (default 1e-10)
% Nmax: maximum # of Chebyshev coefficients
arguments
    f
    I = [-1,1];
    rel_tol = 1e-10;
    Nmax = 100;
end

devalFun = @(x) deval(f, x).';

i = 3;
while true
    test_size = 2^i + 1;
    schbc = chebcoeffs(devalFun, test_size, I);
    amp_k = max(abs(schbc), [], 1);          % 1xN row vector of per-mode peak
    global_max = max(amp_k);
    noise_floor = rel_tol * global_max;
    tail_k = max(abs(schbc(end-4:end, :)), [], 1);

    % Active filter: modes whose peak exceeds noise_floor
    active = amp_k > noise_floor;
    if all(~active)
      relmax = 0;
    else
      relmax = max(tail_k(active) ./ amp_k(active));
    end

    if relmax < rel_tol || (test_size > Nmax)
      thresh = noise_floor;
      [nonzero_row_index, ~] = find(abs(schbc) > thresh);
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
