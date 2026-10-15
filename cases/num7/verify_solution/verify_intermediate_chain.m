function bridge = verify_intermediate_chain(a0,r_initial,abar,params,L,nu,m,h_max,Kmax,terminal,options)
% Sequential full-PDE CAP; no production-cache or proof-file writes.
if nargin<11,options=struct();end
if ~isfield(options,'frame_kind'),options.frame_kind='eigen';end
N=size(a0);if isfield(options,'seed_file'),Cached=load(options.seed_file,'frame');frame=Cached.frame;assert(isequal(frame.m,m)&&sup(norm(frame.V,1))<=1);C=intval(eye(prod(m)))-frame.V*frame.Vi;assert(all(inf(C)<=0 & sup(C)>=0,'all'));else,frame=getting_bridge_frame(abar,m,params,L,nu,options.frame_kind);end;n=prod(m);
state=struct('Phi_acc',intval(eye(n)),'sigma_acc',intval(1), ...
    'gamma',intval(zeros(2,1)),'r_end',intval(r_initial),'a_end',a0(:));
history=struct('h',{},'r',{},'r_end',{},'rho',{},'Pr',{},'Zv',{}, ...
    'schur',{},'neumann',{},'gamma',{},'ba',{});
bridge=struct('verified',false,'steps',0,'duration',intval(0),'rho',intval(Inf), ...
    'state',state,'history',history,'frame',frame,'terminal_test', ...
    test_infinite_radius(Inf,terminal.W,terminal.forcing,terminal.Vq,terminal.anorm), ...
    'reason','No terminal connection has been verified.');
start_step=0;
if isfield(options,'seed_file')
 S=load(options.seed_file);assert(isequal(S.m,m)&&S.h==h_max&&isequal(S.N,N)&&isequal(S.frame,frame)&&isequal(S.state,state));assert(S.test.passed&&sup(S.test.Pr)<0&&sup(S.test.Zv)<1);
 state=S.candidate;history=struct('h',S.h,'r',S.test.radius,'r_end',state.r_end,'rho',getting_bridge_matching(state,N,abar,nu),'Pr',S.test.Pr,'Zv',S.test.Zv,'schur',S.detail.schur,'neumann',S.detail.neumann,'gamma',state.gamma,'ba',S.seg.ba);
 bridge.steps=1;bridge.duration=intval(S.h);bridge.state=state;bridge.history=history;start_step=1;
 stepdir=fullfile(fileparts(mfilename('fullpath')),'middle_steps');if ~isfolder(stepdir),mkdir(stepdir);end
 step=1;duration=bridge.duration;last_step=history;save(fullfile(stepdir,'accepted_step_0001.mat'),'state','duration','last_step','step','history','frame','-v7.3');copyfile(options.seed_file,fullfile(stepdir,'original_probe_step_0001.mat'));
 fprintf('ACTUAL_SEED_REUSED step1 h=%.12g duration=%.12g rho=%.12g r_end=%.12g\n',S.h,sup(duration),sup(history.rho),sup(state.r_end));
end
for k=start_step:Kmax
    rho=getting_bridge_matching(state,N,abar,nu);
    final=test_infinite_radius(rho,terminal.W,terminal.forcing,terminal.Vq,terminal.anorm);
    bridge.rho=rho;bridge.terminal_test=final;
    if final.passed,bridge.verified=true;bridge.reason='';break;end
    if k==Kmax,break;end
    passed=false;h=h_max;attempts=cell(1,7);
    for attempt=1:7
        seg=getting_bridge_segment(state.a_end,N,h,m,params,L,nu);
        try
            local=getting_bridge_bounds(seg,h,params,L,nu,frame);
            [candidate,test,detail]=verify_intermediate_step(state,seg,local,frame,r_initial,nu);
            passed=test.passed;
            attempts{attempt}=struct('h',h,'passed',passed,'test',test,'detail',detail, ...
                'Vq_norm',local.Vq_norm,'anorm',seg.banorm);
        catch err
            if ~ismember(err.identifier,{'getting_bridge_bounds:Schur','getting_int_semi_matrix:Neumann'})
                rethrow(err);
            end
            passed=false;
            attempts{attempt}=struct('h',h,'passed',false,'identifier',err.identifier,'reason',err.message);
        end
        record_root=getenv('GE_VERIFIED_RUN_DIR');assert(~isempty(record_root));stepdir=fullfile(record_root,'middle_steps');if ~isfolder(stepdir),mkdir(stepdir);end
        proof_attempt=struct('prior_state',state,'segment',seg,'result',attempts{attempt},'frame',frame,'duration_before',bridge.duration,'step',k+1,'attempt',attempt);
        if isfield(attempts{attempt},'test'),proof_attempt.local=local;proof_attempt.candidate=candidate;end
        proof_attempt.producer_sha256=ge_file_sha256([mfilename('fullpath') '.m']);proof_attempt.source=ge_source_manifest(getenv('AI_CASE_ROOT'));
        save(fullfile(stepdir,sprintf('step_%04d_attempt_%02d.mat',k+1,attempt)),'proof_attempt','-v7.3');
        if passed,break;end
        h=h/2;
    end
    if ~passed,bridge.reason=sprintf('Local validation failed at intermediate step %d.',k+1);bridge.failed_attempts=attempts;break;end
    state=candidate;bridge.steps=k+1;bridge.duration=bridge.duration+intval(h);
    history(end+1)=struct('h',h,'r',test.radius,'r_end',state.r_end, ...
        'rho',getting_bridge_matching(state,N,abar,nu), ...
        'Pr',test.Pr,'Zv',test.Zv,'schur',detail.schur,'neumann',detail.neumann, ...
        'gamma',state.gamma,'ba',seg.ba); %#ok<AGROW>
    if mod(k+1,100)==0
        fprintf('  intermediate CAP %d: time=%.4f r=%.3e endpoint=%.3e\n', ...
            k+1,sup(bridge.duration),test.radius,sup(state.r_end));
    end
    if isfield(options,'progress_file')
        fprintf('  bridge %d h=%.6g duration=%.6g rho=%.6g r_end=%.6g\n',k+1,h,sup(bridge.duration),sup(history(end).rho),sup(state.r_end));
        duration=bridge.duration;last_step=history(end);step=k+1;
        save(options.progress_file,'state','duration','last_step','step','options','history','attempts','frame','-v7.3');
        save(fullfile(stepdir,sprintf('accepted_step_%04d.mat',step)),'state','duration','last_step','step','history','attempts','frame','-v7.3');
    end
end
bridge.state=state;bridge.history=history;
end
