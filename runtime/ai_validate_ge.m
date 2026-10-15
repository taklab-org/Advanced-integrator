function ai_validate_ge(out)
% Recompute the algebraic step tests from freshly generated local bounds.
R=load(fullfile(out,'result.mat'));F=load(fullfile(out,'result.mat.forward.mat'));
assert(all(structfun(@(x)logical(x),R.validation_status)));
if R.K_bridge>0
 frame=R.bridge.frame;n=prod(frame.m);
 state=struct('Phi_acc',intval(eye(n)),'sigma_acc',intval(1), ...
  'gamma',intval(zeros(2,1)),'r_end',intval(F.re_F),'a_end',F.a_F_mat(:));
 for k=1:R.K_bridge
  files=dir(fullfile(out,'middle_steps',sprintf('step_%04d_attempt_*.mat',k)));
  accepted=false;
  for j=1:numel(files)
   q=load(fullfile(files(j).folder,files(j).name));v=q.proof_attempt;
   assert(isequaln(v.prior_state,state));
   if isfield(v.result,'test')
    [next,test,detail]=verify_intermediate_step(state,v.segment,v.local,frame,intval(F.re_F),F.nu);
    assert(test.passed==v.result.passed);
    if test.passed
     assert(sup(test.Pr)<0&&sup(test.Zv)<1&&detail.schur>0&&detail.neumann<1);
     % Comparison is within one run. Across machines, compare inequalities.
     assert(isequaln(next,v.candidate));accepted=true;
    end
   end
  end
  assert(accepted);state=next;
 end
 rho=getting_bridge_matching(state,size(F.a_F_mat),F.a_star_full,F.nu);
else
 rho=R.bridge.rho;
end
T=test_infinite_radius(rho,F.terminal.W,F.terminal.forcing,F.terminal.Vq,F.terminal.anorm);
assert(T.passed&&sup(T.Pr)<0&&sup(T.Zv)<1);
ai_write_json(fullfile(out,'validation.json'),struct('passed',true,'steps',R.K_bridge, ...
 'terminal_Pr_upper',sup(T.Pr),'terminal_Z_upper',sup(T.Zv)));
end
