function Un_endpt = getting_Un_endpt(Phi,r1)
%
Phi(2:end,:,:) = 2*Phi(2:end,:,:); % convert to two-sided chebyshev
% permute(sum(intval(Phi),1),[2,3,1])
Un_endpt = (permute(sum(intval(Phi),1),[2,3,1]) + r1*infsup(-1,1));
