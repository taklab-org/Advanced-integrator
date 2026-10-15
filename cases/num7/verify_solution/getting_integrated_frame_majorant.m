function [bound,evidence] = getting_integrated_frame_majorant(data,h,left,right)
% Bound the inhomogeneous Volterra operator in a fixed frame.
% J integrates |left*E(t)*right| before discarding signed cancellation.
% K remains the original checked Duhamel defect bound.
n=size(data.V,1);assert(inf(intval(h))>=0 && sup(intval(h))<=sup(data.h));assert(sup(norm(data.K,inf))<1);
[~,~,J0]=getting_matrix_exp_bound(data,h);
[~,~,Jframe]=getting_matrix_exp_bound(data,h,left,right);
resolved=data.Bneu\J0;
correction=abs(intval(left))*data.K*resolved*abs(intval(right));
bound=intval(max(0,sup(Jframe+correction)));
evidence=struct('J0',J0,'Jframe',Jframe,'resolved',resolved,'correction',correction,'K',data.K,'Bneu',data.Bneu);
end
