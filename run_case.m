function run_case(name,intlabRoot,out,options)
%RUN_CASE Recompute one published example without any precomputed MAT files.
% run_case('num7',getenv('INTLAB_ROOT'),'/absolute/output',struct('Workers',4))
% Use a fresh MATLAB process; never run two examples in one MATLAB session.
if nargin<4,options=struct;end
root=fileparts(mfilename('fullpath'));
assert(ismember(name,{'num5','num7','pbvp'}),'Unknown example.');
if nargin<2||isempty(intlabRoot),intlabRoot=getenv('INTLAB_ROOT');end
if nargin<3||isempty(out),out=fullfile(root,'outputs',[name '_' datestr(now,'yyyymmdd_HHMMSS')]);end
if ~isfolder(out),mkdir(out);end
out=char(java.io.File(out).getCanonicalPath());
entries=dir(out);assert(isempty(setdiff({entries.name},{'.','..'})), ...
 'Choose a new, empty output directory; no precomputed data are accepted.');
assert(~isfile(fullfile(out,'RUN_STARTED')),'Choose a new output directory.');
options=defaults(options);
assert(options.Workers>=1&&options.Workers<=10&&options.Threads>=1&&options.ClientThreads>=1);
addpath(fullfile(root,'runtime'),root,'-begin');
oldpath=getenv('MATLABPATH');oldpwd=pwd;
cleanup=onCleanup(@()finish(oldpath,oldpwd));
setenv('MATLABPATH',[fullfile(root,'runtime') pathsep root]);
setenv('GE_VERIFIED_RUN_DIR',out);setenv('AI_CASE_ROOT',fullfile(root,'cases',name));
setenv('AI_OPTIONS',jsonencode(options));setenv('AI_SMOKE',num2str(options.Smoke));
write_text(fullfile(out,'RUN_STARTED'),datestr(now,30));
diary(fullfile(out,'run.log'));elapsed=tic;
try
 sourceManifest=ai_release_manifest(root);
 ai_write_json(fullfile(out,'source-manifest.json'),sourceManifest);
 maxNumCompThreads(options.ClientThreads);
 global INTLAB_CONST
 if isempty(INTLAB_CONST)
  restoredefaultpath;addpath(fullfile(root,'runtime'),root,'-begin');
  clientReport=intlab_client_start(intlabRoot);
 else
  clientReport=intlab_runtime_check();
  assert(strcmp(char(java.io.File(clientReport.intlabRoot).getCanonicalPath()), ...
    char(java.io.File(intlabRoot).getCanonicalPath())),'An unrelated INTLAB is active.');
 end
 clientChecks=ge_check();
 save(fullfile(out,'client_environment.mat'),'clientReport','clientChecks');
 ai_add_case_paths(getenv('AI_CASE_ROOT'));addpath(fullfile(root,'runtime'),'-begin');
 save(fullfile(out,'configuration.mat'),'name','options','clientReport');
 ai_write_json(fullfile(out,'configuration.json'),struct('example',name,'options',options, ...
  'MATLAB',version,'BLAS',version('-blas'),'platform',computer));
 if strcmp(name,'pbvp')
  % PBVP retains its recorded Chebyshev kernel; the independent-column
  % shadow belongs to the GE cases only.
  addpath(fullfile(getenv('AI_CASE_ROOT'),'defect','Chebyshev'),'-begin');
  ai_generate_pbvp(out);
  cd(out);script_proof_connexion4to3;
  assert(bvp_ok&&all(sup(Y_bvp+Zr*max(intval(r)))<r)&&max(sup(Zr))<1);
  save(fullfile(out,'result.mat'),'r','bvp_ok','Y_bvp','Zr','M','G_mat','Z0_mat', ...
   'Z1_mat','B_mat','U_DP','err_lin','g_theta','g_phi','g_u_tau','g_u_X','P4', ...
   'tau_vec','n_cheb_vec','K','r1_vec','theta','phi','-v7.3');
  summary=struct('success',bvp_ok,'steps',K,'radii',r(:).','Z_upper',sup(Zr(:)).', ...
   'margins',inf(intval(r)-Y_bvp-Zr*max(intval(r))).');
 else
  if strcmp(name,'num5'),S=initial_num5();N=20;else,S=initial_num7();N=28;end
  ic=fullfile(out,'initial.mat');save(ic,'-struct','S');
  proof_GE_2DOK_core(ic,N,fullfile(out,'phase2.mat'),fullfile(out,'phase3.mat'), ...
   fullfile(out,'result.mat'),name);
  R=load(fullfile(out,'result.mat'));
  if options.Smoke
   assert(R.source_data_revalidated);summary=struct('success',true,'smoke_only',true);
  else
   assert(R.success&&all(structfun(@(x)logical(x),R.validation_status))&&sup(R.Pr)<0&&sup(R.Zv)<1);
   if strcmp(name,'num7'),assert(R.K_F==128&&R.tau_F==8);else,assert(R.K_F==208&&R.tau_F==13);end
   ai_validate_ge(out);
   summary=struct('success',R.success,'validation_status',R.validation_status,'tau_F',R.tau_F, ...
    'steps',R.K_F,'bridge_steps',R.K_bridge,'entry_time',R.tau_total, ...
    'terminal_radius',R.r_up,'Pr_upper',sup(R.Pr),'Z_upper',sup(R.Zv));
  end
 end
 assert(isequal(sourceManifest,ai_release_manifest(root)),'Release sources changed during calculation.');
 summary.elapsed_seconds=toc(elapsed);ai_write_json(fullfile(out,'summary.json'),summary);
 write_text(fullfile(out,'SUCCESS'),'Strict inequalities checked; see summary.json.');
catch ME
 write_text(fullfile(out,'FAILED'),getReport(ME,'extended','hyperlinks','off'));
 rethrow(ME);
end
diary off;
end
function o=defaults(o)
d=struct('Workers',2,'Threads',2,'ClientThreads',2,'Smoke',false);
f=fieldnames(o);assert(all(ismember(f,fieldnames(d))),'Unknown option.');
for f=fieldnames(d).',if ~isfield(o,f{1}),o.(f{1})=d.(f{1});end,end
end
function finish(p,d)
if exist('gcp','file')&&~isempty(gcp('nocreate')),delete(gcp('nocreate'));end
setenv('MATLABPATH',p);cd(d);diary off;
end
function write_text(f,s)
fid=fopen(f,'w');assert(fid>0);c=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',s);
end
