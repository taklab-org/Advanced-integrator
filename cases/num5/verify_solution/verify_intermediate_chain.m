function bridge = verify_intermediate_chain(a0,r_initial,abar,params,L,nu,m,h_max,Kmax,terminal,options)
% Sequential full-PDE CAP; no production-cache or proof-file writes.
if nargin<11,options=struct();end
if ~isfield(options,'frame_kind'),options.frame_kind='eigen';end
N=size(a0);frame=getting_bridge_frame(abar,m,params,L,nu,options.frame_kind);n=prod(m);
state=struct('Phi_acc',intval(eye(n)),'sigma_acc',intval(1), ...
    'gamma',intval(zeros(2,1)),'r_end',intval(r_initial),'a_end',a0(:));
history=struct('h',{},'r',{},'r_end',{},'rho',{},'Pr',{},'Zv',{}, ...
    'schur',{},'neumann',{},'gamma',{},'ba',{});
bridge=struct('verified',false,'steps',0,'duration',intval(0),'rho',intval(Inf), ...
    'state',state,'history',history,'frame',frame,'terminal_test', ...
    test_infinite_radius(Inf,terminal.W,terminal.forcing,terminal.Vq,terminal.anorm), ...
    'reason','No terminal connection has been verified.');
for k=0:Kmax
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
        save(options.progress_file,'state','duration','last_step','step','options','-v7.3');
    end
end
bridge.state=state;bridge.history=history;
end
