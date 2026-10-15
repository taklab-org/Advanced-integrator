function report = intlab_worker_start(clientReport)
%INTLAB_WORKER_START Minimum initialization against a prebuilt client cache.
% For fresh process workers with automatic full startup suppressed beforehand.
% clientReport must come from intlab_runtime_check on the same runtime.
assert(isempty(getenv('FLML_TRUSTED_BOOTSTRAP_RECORD')), ...
    'intlabStartup:Formal','Use the formal project worker startup owner.');
assert(strcmp(version,clientReport.matlabVersion) ...
    && strcmp(computer,clientReport.computer), ...
    'intlabStartup:Platform','Prepare the cache with this MATLAB version and platform.');
global INTLAB_CONST %#ok<GVMIS>
assert(isempty(INTLAB_CONST), ...
    'intlabStartup:WorkerState', ...
    'Worker startup already initialized INTLAB; inspect startup.m before creating a fresh pool.');
runtimeRoot = clientReport.intlabRoot;
minimumPath = fullfile(runtimeRoot,'startintlab_minimum.m');
assert(isfile(minimumPath) && isfile(clientReport.cachePath), ...
    'intlabStartup:MissingCache','A reviewed minimum initializer and prebuilt cache are required.');
cacheBefore = readCache(clientReport.cachePath);
cacheMetaBefore = dir(clientReport.cachePath);
saved = load(clientReport.cachePath,'INTLAB_CONST');
assert(isfield(saved,'INTLAB_CONST') && isstruct(saved.INTLAB_CONST) ...
    && isfield(saved.INTLAB_CONST,'ROUNDING_SUCCESS') ...
    && isscalar(saved.INTLAB_CONST.ROUNDING_SUCCESS) ...
    && saved.INTLAB_CONST.ROUNDING_SUCCESS == 1 ...
    && isfield(saved.INTLAB_CONST,'INTLABPATH') ...
    && sameDirectory(saved.INTLAB_CONST.INTLABPATH,runtimeRoot), ...
    'intlabStartup:Cache','The cache is invalid or belongs to a different runtime root.');
initialDirectory = pwd;
restoreDirectory = onCleanup(@() cd(initialDirectory));
addpath(runtimeRoot,'-begin');
assert(strcmp(which('startintlab'),fullfile(runtimeRoot,'startintlab.m')) ...
    && strcmp(which('startintlab_minimum'),minimumPath), ...
    'intlabStartup:Shadow','An INTLAB startup entry is shadowed.');
invokeMinimumStartup();
report = intlab_runtime_check();
cacheMetaAfter = dir(clientReport.cachePath);
assert(strcmp(report.cachePath,clientReport.cachePath) ...
    && isequal(readCache(clientReport.cachePath),cacheBefore) ...
    && cacheMetaAfter.bytes == cacheMetaBefore.bytes ...
    && cacheMetaAfter.datenum == cacheMetaBefore.datenum, ...
    'intlabStartup:CacheChanged','The prepared cache changed during worker initialization.');
report.cacheUnchanged = true;
end

function invokeMinimumStartup()
% Isolate caller-workspace variables if the initializer is a local script.
startintlab_minimum;
end

function bytes = readCache(cachePath)
fid = fopen(cachePath,'rb');
assert(fid ~= -1,'intlabStartup:CacheRead','Cannot read the prepared INTLAB cache.');
closeFile = onCleanup(@() fclose(fid));
bytes = fread(fid,Inf,'*uint8');
assert(~isempty(bytes),'intlabStartup:EmptyCache','The prepared cache is empty.');
end

function same = sameDirectory(a,b)
% JVM canonicalization recognizes symlink aliases (including Dropbox paths).
if usejava('jvm')
    a = char(java.io.File(a).getCanonicalPath());
    b = char(java.io.File(b).getCanonicalPath());
else
    [okA,infoA] = fileattrib(a); [okB,infoB] = fileattrib(b);
    assert(okA && okB,'intlabStartup:Path','Cannot resolve INTLAB directories.');
    a = infoA.Name; b = infoB.Name;
end
if ispc, same = strcmpi(a,b); else, same = strcmp(a,b); end
end
