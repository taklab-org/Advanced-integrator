function cc = pdecheb2_relative_v3(f, I, rel_tol, abs_floor, Nmax)
% pdecheb2 with HYBRID noise filter:
%   noise_filter = max(rel_tol * global_max_amp, abs_floor)
% Active mode: amp_k > noise_filter.
% Stop: max(tail_k(active) ./ amp_k(active)) < rel_tol.
arguments
    f
    I = [-1,1];
    rel_tol = 1e-14;
    abs_floor = 1e-13;
    Nmax = 100;
end

devalFun = @(x) deval(f, x).';
i = 3;
while true
    test_size = 2^i + 1;
    schbc = chebcoeffs(devalFun, test_size, I);
    amp_k = max(abs(schbc), [], 1);
    global_max = max(amp_k);
    noise_filter = max(rel_tol * global_max, abs_floor);
    tail_k = max(abs(schbc(end-4:end, :)), [], 1);

    active = amp_k > noise_filter;
    if all(~active)
      relmax = 0;
    else
      relmax = max(tail_k(active) ./ amp_k(active));
    end
    if relmax < rel_tol || (test_size > Nmax)
      [nonzero_row_index, ~] = find(abs(schbc) > noise_filter);
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
