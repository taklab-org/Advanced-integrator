function summaries = reproduce(intlabRoot,outputRoot,example,options)
%REPRODUCE Run examples sequentially, each in a fresh MATLAB process.
% reproduce('/path/to/INTLAB','/absolute/new/output')
% reproduce('/path/to/INTLAB','/absolute/new/output','num7', ...
%     struct('Workers',4,'Threads',2,'ClientThreads',8))
% options.DryRun=true checks/plans launches without starting calculations.
if nargin<1||isempty(intlabRoot),intlabRoot=getenv('INTLAB_ROOT');end
if nargin<2||isempty(outputRoot)
 outputRoot=fullfile(fileparts(mfilename('fullpath')),'outputs',datestr(now,'yyyymmdd_HHMMSS'));
end
if nargin<3||isempty(example),example='all';end
if nargin<4,options=struct;end
root=fileparts(mfilename('fullpath'));
intlabRoot=char(intlabRoot);example=char(example);
assert(~isempty(intlabRoot)&&isfolder(intlabRoot),'AdvancedIntegrator:INTLAB', ...
 'Pass the installed INTLAB directory, or set INTLAB_ROOT.');
assert(ismember(example,{'all','num5','num7','pbvp'}),'AdvancedIntegrator:Example','Unknown example.');
if ispc,exe=fullfile(matlabroot,'bin','matlab.exe');else,exe=fullfile(matlabroot,'bin','matlab');end
d=struct('Workers',2,'Threads',2,'ClientThreads',2,'Smoke',false, ...
 'MatlabExecutable',exe,'DryRun',false);
assert(isstruct(options)&&isscalar(options),'Options must be a scalar struct.');
assert(all(ismember(fieldnames(options),fieldnames(d))),'Unknown option.');
for f=fieldnames(d).',if ~isfield(options,f{1}),options.(f{1})=d.(f{1});end,end
for f={'Workers','Threads','ClientThreads'}
 validateattributes(options.(f{1}),{'numeric'},{'scalar','integer','positive','finite'});
end
assert(options.Workers<=10,'At most ten workers are supported.');
validateattributes(options.Smoke,{'logical','numeric'},{'scalar','binary'});
validateattributes(options.DryRun,{'logical','numeric'},{'scalar','binary'});
assert(isfile(options.MatlabExecutable),'MATLAB executable not found.');
outputRoot=char(java.io.File(char(outputRoot)).getCanonicalPath());
intlabRoot=char(java.io.File(intlabRoot).getCanonicalPath());
if strcmp(example,'all'),names={'num5','num7','pbvp'};else,names={example};end
runOptions=rmfield(options,{'MatlabExecutable','DryRun'});
% Preflight every output before launching the first long calculation.
for k=1:numel(names)
 out=fullfile(outputRoot,names{k});
 if isfolder(out)
  f=dir(out);assert(isempty(setdiff({f.name},{'.','..'})), ...
   'AdvancedIntegrator:OutputNotEmpty','Output is not empty: %s',out);
 else
  assert(~isfile(out),'Output path is a file: %s',out);
 end
end
summaries=cell(1,numel(names));
for k=1:numel(names)
 name=names{k};out=fullfile(outputRoot,name);
 fprintf('Advanced-integrator: %s -> %s\n',name,out);
 if options.DryRun
  summaries{k}=struct('example',name,'output',out,'options',runOptions,'dry_run',true);
  continue
 end
 if ~isfolder(outputRoot),mkdir(outputRoot);end
 if ~isfolder(out),mkdir(out);end
 % Java ProcessBuilder passes argv directly: no shell quoting of user paths,
 % no permanent changes to the launching MATLAB session's environment.
 expression=['addpath(getenv(''AI_REPOSITORY''));' ...
  'run_case(getenv(''AI_EXAMPLE''),getenv(''INTLAB_ROOT''),' ...
  'getenv(''AI_OUTPUT''),jsondecode(getenv(''AI_LAUNCH_OPTIONS'')))'];
 args=javaArray('java.lang.String',3);
 args(1)=java.lang.String(char(options.MatlabExecutable));
 args(2)=java.lang.String('-batch');args(3)=java.lang.String(expression);
 builder=java.lang.ProcessBuilder(args);
 builder.directory(java.io.File(root));
 env=builder.environment();
 env.put('INTLAB_ROOT',intlabRoot);env.put('AI_REPOSITORY',root);
 env.put('AI_EXAMPLE',name);env.put('AI_OUTPUT',out);
 env.put('AI_LAUNCH_OPTIONS',jsonencode(runOptions));
 env.put('MATLABPATH',[fullfile(root,'runtime') pathsep root]);
 prefs=fullfile(outputRoot,[name '_prefs']);if ~isfolder(prefs),mkdir(prefs);end
 env.put('MATLAB_PREFDIR',prefs);
 log=fullfile(outputRoot,[name '_launcher.log']);
 builder.redirectErrorStream(true);builder.redirectOutput(java.io.File(log));
 process=builder.start();
 fprintf('  Fresh MATLAB started. Log: %s\n',log);
 % A blocking numerical run; never launch the next example before exit.
 while process.isAlive(),pause(.25);end
 status=process.exitValue();
 assert(status==0,'AdvancedIntegrator:ChildFailed', ...
  '%s MATLAB exited with status %d. See %s',name,status,log);
 assert(isfile(fullfile(out,'SUCCESS'))&&isfile(fullfile(out,'summary.json')), ...
  'AdvancedIntegrator:MissingResult','No successful result for %s. See %s',name,log);
 summary=jsondecode(fileread(fullfile(out,'summary.json')));
 assert(summary.success,'AdvancedIntegrator:FailedResult','Result did not pass: %s',name);
 summaries{k}=summary;disp(summary);
end
end
