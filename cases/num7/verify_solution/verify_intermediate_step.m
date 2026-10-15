function [next,test,detail] = verify_intermediate_step(state,seg,local,frame,r_initial,nu)
% Initial error stays in one JOINT physical finite/tail 1-ball.
% Phi_acc retains the finite correlations. gamma has only TWO components:
% a finite-frame 1-ball radius and a tail 1-norm radius.
n=prod(frame.m);next=state;
P=intval(zeros(n+1));P(1:n,1:n)=abs(frame.Vi*state.Phi_acc);P(end,end)=state.sigma_acc;
[Sphys,~]=getting_finite_flow(local.finite,local.h,frame.D,frame.Di*state.Phi_acc);
[Sz,~]=getting_finite_flow(local.finite,local.h,frame.Vi*frame.D,frame.Di*state.Phi_acc);
bare_phys=intval(zeros(n+1));bare_phys(1:n,1:n)=Sphys;bare_phys(end,end)=state.sigma_acc;
bare_z=intval(zeros(n+1));bare_z(1:n,1:n)=Sz;bare_z(end,end)=state.sigma_acc;
coupling_initial=local.G*local.B*bare_z;
known=norm(bare_phys,1)+norm(coupling_initial,1);
if isfield(local,'tail_initial_correction')
 assert(sup(state.sigma_acc)<=1);
 tail_correction=local.tail_initial_correction*state.sigma_acc;
 refined=norm(bare_phys,1)+max(norm(coupling_initial(:,1:n),1),norm(tail_correction,1));
 known=intval(min(sup(known),sup(refined)));
end
gap0=getting_bridge_gap(seg.start_error,frame,nu);
Y_initial=intval(r_initial)*known;Y_gamma=sum(local.S2*state.gamma);
Y_gap=sum(local.S*gap0);Y_forcing=sum(local.forcing);
Y=Y_initial+Y_gamma+Y_gap+Y_forcing;
test=test_cubic_radius(Y,local.Vq_norm,seg.banorm);
detail=struct('Y',Y,'known_initial',known,'schur',local.schur,'neumann',local.neumann);
detail.Y_initial=Y_initial;detail.Y_gamma=Y_gamma;detail.Y_gap=Y_gap;detail.Y_forcing=Y_forcing;
if ~test.passed,return;end
r=intval(test.radius);Nterm=(3*seg.banorm+r)*r^2;
Rknown=getting_block_norms(local.R*P,n);
known_remainder=intval(r_initial)*max(Rknown,[],2);
if isfield(local,'tail_initial_correction')
 original_columns=local.R*P;tail_columns=local.tail_initial_correction*state.sigma_acc;
 finite_groups=[max(sum(abs(original_columns(1:n,1:n)),1));max(abs(original_columns(end,1:n)))];
 tail_groups=[max(sum(abs(tail_columns(1:n,:)),1));max(abs(tail_columns(end,:)))];
 refined_remainder=intval(r_initial)*max(finite_groups,tail_groups);
 known_remainder=intval(min(sup(known_remainder),sup(refined_remainder)));
end
gap1=getting_bridge_gap(seg.end_error,frame,nu);
gap_end=local.H2*[sum(gap0(1:n));gap0(end)];
gamma=local.H2*state.gamma+known_remainder+gap_end ...
    +[sum(local.forcing(1:n));local.forcing(end)] ...
    +max(local.Vq2,[],2)*Nterm+[sum(gap1(1:n));gap1(end)];
next.gamma=intval(max(0,sup(gamma)));
next.Phi_acc=local.Phi*state.Phi_acc;
next.sigma_acc=intval(sup(local.sigma*state.sigma_acc));
next.r_end=intval(r_initial)*max(norm(next.Phi_acc,1),next.sigma_acc)+sum(next.gamma);
% The local ball gives a second, independent endpoint bound.
next.r_end=intval(min(sup(next.r_end),sup(r+wnorm(seg.end_error,intval(nu)))));
next.a_end=seg.a_end;
detail.Pr=test.Pr;detail.Zv=test.Zv;detail.gamma=next.gamma;detail.r_end=next.r_end;
end
