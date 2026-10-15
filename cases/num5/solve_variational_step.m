function [bUn, bUn_endpt, r1, verification] = solve_variational_step(ba, h, m, params, L)
% Wrapper for parfeval: solve one variational step for OK and return results
  [bUn, r1, verification] = verify_variationaleq(ba, h, m, params, L);
  bUn_endpt = getting_Un_endpt(bUn, r1);
  assert(all(isfinite(sup(abs(bUn_endpt))),'all'),'GE:VariationalEndpoint','Nonfinite endpoint enclosure.');
end
