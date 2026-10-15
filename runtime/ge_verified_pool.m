function reports=ge_verified_pool(varargin)
% Explicit process pool; requested historical counts are replaced by options.
assert(isempty(gcp('nocreate')),'Only this run may own the process pool.');
out=getenv('GE_VERIFIED_RUN_DIR');o=jsondecode(getenv('AI_OPTIONS'));
E=load(fullfile(out,'client_environment.mat'),'clientReport');
c=parcluster('Processes');c.NumWorkers=max(c.NumWorkers,o.Workers);
c.JobStorageLocation=fullfile(out,'jobs');if ~isfolder(c.JobStorageLocation),mkdir(c.JobStorageLocation);end
p=parpool(c,o.Workers,'IdleTimeout',Inf);
caseRoot=getenv('AI_CASE_ROOT');runtime=fileparts(mfilename('fullpath'));
% Full startup fallback is serialized, so workers never write a shared cache
% concurrently. A reviewed minimum initializer is used if installed.
spmd
 restoredefaultpath;addpath(runtime,'-begin');
 ai_add_case_paths(caseRoot);addpath(runtime,'-begin');
 for owner=1:numlabs
  if labindex==owner,record=ai_worker_init(E.clientReport,o.Threads);end
  labBarrier;
 end
end
reports=cell(1,p.NumWorkers);for k=1:p.NumWorkers,reports{k}=record{k};end
assert(numel(reports)==p.NumWorkers);
save(fullfile(out,'workers_verified.mat'),'reports');
save(fullfile(out,['workers_' datestr(now,'yyyymmdd_HHMMSS') '.mat']),'reports');
end
