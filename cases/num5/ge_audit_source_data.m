function [passed,record]=ge_audit_source_data(ic_file,N_requested,p2_file,p3_file,root)
passed=false;record=struct('passed',false,'reason','Legacy files have no checked generation record.');
w2=whos('-file',p2_file);w3=whos('-file',p3_file);
if ~ismember('generation_record',{w2.name}) || ~ismember('generation_record',{w3.name}),return;end
P=load(p2_file);Q=load(p3_file);S=load(ic_file);g2=P.generation_record;g3=Q.generation_record;
assert(g2.schema==1 && g3.schema==1 && g2.complete && g3.complete && g2.phase==2 && g3.phase==3);
source=ge_source_manifest(root);
assert(isequal(source,g2.source) && isequal(source,g3.source),'GE:SourceMismatch','Generated data and executing sources differ.');
assert(strcmp(g2.ic_sha256,ge_file_sha256(ic_file)),'GE:InitialMismatch','Initial condition changed.');
p2_hash=ge_file_sha256(p2_file);assert(strcmp(g3.phase2_sha256,p2_hash),'GE:PhaseMismatch','Phase 3 belongs to different Phase 2 data.');
assert(isequal(P.N,[N_requested,N_requested]) && isequal(P.N,g2.N));
assert(isequal(P.params,S.params) && isequal(P.L,S.L));
ni=round(sqrt(numel(S.a0_full)));initial=reshape(S.a0_full,ni,ni);padded=zeros(P.N);o=min(ni,N_requested);
padded(1:o,1:o)=initial(1:o,1:o);initial(1:o,1:o)=0;
assert(all(initial==0,'all'),'GE:InitialTruncation','The requested N omits nonzero initial modes.');
assert(isequal(P.a0_init_vec{1}(:),padded(:)));
K=P.num_integration;assert(K==Q.num_integration && K==g2.steps && K==g3.steps);
assert(numel(P.ba_vec)==K && numel(P.a0_init_vec)==K+1 && numel(Q.bUn_vec)==K && numel(Q.bUn_endpt_vec)==K);
assert(all(P.tau_vec>0) && sum(P.tau_vec)==g2.target_time && all(g2.schur_margins>0));
assert(all(isfinite(sup(Q.r1_vec)) & inf(Q.r1_vec)>=0,'all'));
margin=Inf;
for k=1:K
 v=g3.verification{k};assert(v.passed && all(inf(v.margin)>0,'all') && isequal(v.m,P.m) && v.h==P.tau_vec(k));
 margin=min(margin,min(inf(v.margin),[],'all'));
 assert(all(isfinite(sup(abs(Q.bUn_endpt_vec{k}))),'all'));
end
run=getenv('GE_VERIFIED_RUN_DIR');assert(~isempty(run),'GE:RuntimeEvidence','Verified runtime required.');
E=load(fullfile(run,'client_environment.mat'),'clientChecks');
assert(E.clientChecks.runtime.smokePassed && E.clientChecks.intervalMatrixPassed && E.clientChecks.pid==feature('getpid'));
W=load(fullfile(run,'workers_verified.mat'),'reports');reports=W.reports;
assert(~isempty(reports));pids=zeros(numel(reports),1);
for k=1:numel(reports)
 if iscell(reports),v=reports{k};else,v=reports(k);end
 assert(v.startup.smokePassed && v.startup.cacheUnchanged && v.check.intervalMatrixPassed);
 pids(k)=v.check.pid;
end
assert(numel(unique(pids))==numel(pids));
passed=true;record=struct('passed',true,'reason','Fresh source-bound Phase 2/3 data and all generation checks passed.', ...
 'ic_sha256',g2.ic_sha256,'phase2_sha256',p2_hash,'phase3_sha256',ge_file_sha256(p3_file), ...
 'source',source,'N',P.N,'steps',K,'min_variational_margin',margin,'worker_count',numel(reports), ...
 'client_checks',E.clientChecks);
fprintf('SOURCE_DATA_AUDIT PASSED: N=[%d,%d], %d segments, min variational margin %.6g\n',P.N,K,margin);
end
