function [S,Phi] = getting_finite_flow(data,h,left,right)
% Finite non-autonomous flow inclusion from the checked Neumann defect.
n=size(data.V,1);if nargin<3,left=intval(eye(n));end
if nargin<4,right=intval(eye(n));end
if isfield(data,'h'),assert(sup(intval(h))<=sup(data.h),'Flow bound exceeds its verified window.');end
if ~isfield(data,'error')
    Sbar=getting_matrix_exp_bound(data,h);
    Strue=data.Bneu\Sbar;
    data.error=data.K*Strue;
end
[Sbar,center]=getting_matrix_exp_bound(data,h,left,right);
radius=abs(intval(left))*data.error*abs(intval(right));
S=Sbar+radius;
% The OK variational flow is real, even if its approximate diagonalization
% uses complex coordinates. Project the center to its real part.
Phi=real(center)+midrad(zeros(size(center)),sup(radius));
end
