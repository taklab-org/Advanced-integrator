function resolved = getting_time_forcing(f,m,h,params,L,nu,B,Vfin_w)
% A time-dependent forcing uses its per-mode Chebyshev coefficient majorant.
% Reuse the ordinary, mode-resolved Volterra bound, not a new inverse bound.
f=intval(f);mode_bound=reshape(abs(f(1,:,:))+2*sum(abs(f(2:end,:,:)),1),size(f,2),size(f,3));
resolved=getting_constant_forcing(mode_bound,m,h,params,L,nu,B,Vfin_w);
end
