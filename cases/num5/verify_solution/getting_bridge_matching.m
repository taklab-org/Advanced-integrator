function rho = getting_bridge_matching(state,N,abar,nu)
% Full-support distance from the nominal endpoint to the terminal center.
common=max(N,size(abar));e=intval(zeros(common));
e(1:N(1),1:N(2))=reshape(intval(state.a_end),N);
e(1:size(abar,1),1:size(abar,2))=e(1:size(abar,1),1:size(abar,2))-intval(abar);
rho=state.r_end+wnorm(e,intval(nu));
end
