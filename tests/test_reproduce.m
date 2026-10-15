function test_reproduce()
%TEST_REPRODUCE Check planning, isolation and real child failure propagation.
% Does not start numerical proof work or initialize INTLAB.
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
work=tempname;mkdir(work);cleanup=onCleanup(@()note_retained(work));
fakeIntlab=fullfile(work,'INTLAB directory with spaces');mkdir(fakeIntlab);
out=fullfile(work,'output with spaces');oldpath=path;oldenv=getenv('MATLABPATH');
oldroot=getenv('INTLAB_ROOT');oldpwd=pwd;
s=reproduce(fakeIntlab,out,'all',struct('DryRun',true));
assert(numel(s)==3&&strcmp(s{2}.example,'num7')&&all(cellfun(@(x)x.dry_run,s)));
assert(~isfolder(out));assert_environment(oldpath,oldenv,oldroot,oldpwd);
% Validate every destination up front; num5 must not start if num7 is occupied.
mkdir(out);mkdir(fullfile(out,'num7'));fid=fopen(fullfile(out,'num7','existing'),'w');fclose(fid);
expect_error(@()reproduce(fakeIntlab,out),'AdvancedIntegrator:OutputNotEmpty');
assert(~isfolder(fullfile(out,'num5')));
expect_error(@()reproduce(fakeIntlab,fullfile(work,'bad'),'unknown'),'AdvancedIntegrator:Example');
% Launch an actual fresh MATLAB. An empty directory is deliberately not an
% INTLAB installation: the child's nonzero exit must never become success.
childout=fullfile(work,'child output with spaces');
expect_error(@()reproduce(fakeIntlab,childout,'pbvp'),'AdvancedIntegrator:ChildFailed');
assert(isfile(fullfile(childout,'pbvp_launcher.log')));
assert(isfile(fullfile(childout,'pbvp','FAILED')));
assert(~isfile(fullfile(childout,'pbvp','SUCCESS')));
log=fileread(fullfile(childout,'pbvp_launcher.log'));
assert(contains(log,'startintlab')||contains(log,'intlabStartup:Root'));
assert_environment(oldpath,oldenv,oldroot,oldpwd);
% Positive sequence test: a private fixture replaces only the numerical
% worker entry. These SUCCESS files certify the test fixture, not a proof.
fixture=fullfile(work,'fixture repository');mkdir(fixture);mkdir(fullfile(fixture,'runtime'));
copyfile(fullfile(root,'reproduce.m'),fixture);
copyfile(fullfile(root,'runtime','startup.m'),fullfile(fixture,'runtime'));
fid=fopen(fullfile(fixture,'run_case.m'),'w');
fprintf(fid,'%s\n', ...
 'function run_case(name,intlabRoot,out,options)', ...
 'assert(options.Workers==3);', ...
 's=struct(''success'',true,''fixture'',true,''example'',name,''pid'',feature(''getpid''));', ...
 'fid=fopen(fullfile(out,''summary.json''),''w'');fprintf(fid,''%s'',jsonencode(s));fclose(fid);', ...
 'fid=fopen(fullfile(out,''SUCCESS''),''w'');fprintf(fid,''Test fixture only'');fclose(fid);', ...
 'end');fclose(fid);
restore=onCleanup(@()restore_launcher_path(oldpath,oldpwd));
cd(fixture);clear reproduce;addpath(fixture,'-begin');
assert(strcmp(which('reproduce'),fullfile(fixture,'reproduce.m')));
s=reproduce(fakeIntlab,fullfile(work,'fixture output'),'all',struct('Workers',3));
assert(numel(s)==3&&all(cellfun(@(x)x.success&&x.fixture,s)));
ids=cellfun(@(x)x.pid,s);assert(numel(unique(ids))==3&&all(ids~=feature('getpid')));
assert(isequal(cellfun(@(x)x.example,s,'UniformOutput',false),{'num5','num7','pbvp'}));
clear restore;
assert_environment(oldpath,oldenv,oldroot,oldpwd);
fprintf('MATLAB_LAUNCHER_TESTS_PASSED\n');
remove_test_tree(work);clear cleanup;
end
function restore_launcher_path(p,d)
clear reproduce;path(p);cd(d);
end
function expect_error(f,id)
try,f();catch e,assert(strcmp(e.identifier,id),'Unexpected error: %s',e.message);return;end
error('Expected error %s was not raised.',id);
end
function assert_environment(p,e,r,d)
assert(strcmp(path,p)&&strcmp(getenv('MATLABPATH'),e) ...
 &&strcmp(getenv('INTLAB_ROOT'),r)&&strcmp(pwd,d),'Parent session was modified.');
end
function remove_test_tree(p)
% Only this test's explicit tempname directory is removed after children exit.
if isfolder(p),rmdir(p,'s');end
end
function note_retained(p)
if isfolder(p),fprintf('Failed test artifacts retained: %s\n',p);end
end
