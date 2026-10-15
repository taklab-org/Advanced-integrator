function report = intlab_client_start(intlabRoot)
%INTLAB_CLIENT_START Full initialization for an ordinary serial MATLAB client.
% May create/update INTLAB's cache. Use a dedicated prepared runtime if shared.
% Formal proof launchers own their startup; do not call this before them.
assert(isempty(getenv('FLML_TRUSTED_BOOTSTRAP_RECORD')), ...
    'intlabStartup:Formal','Use the formal project startup owner.');
if exist('getCurrentTask','file')
    assert(isempty(getCurrentTask()),'intlabStartup:WorkerFullStartup', ...
        'Full startup belongs to the serial client; use the prepared-cache worker owner.');
end
if nargin == 0 || isempty(intlabRoot)
    intlabRoot = getenv('INTLAB_ROOT');
end
assert(~isempty(intlabRoot) && isfile(fullfile(intlabRoot,'startintlab.m')), ...
    'intlabStartup:Root','Pass the installed INTLAB root or set INTLAB_ROOT.');
global INTLAB_CONST %#ok<GVMIS>
assert(isempty(INTLAB_CONST), ...
    'intlabStartup:ExistingState', ...
    'INTLAB is already present or partly initialized; check it or use a fresh process.');
initialDirectory = pwd;
restoreDirectory = onCleanup(@() cd(initialDirectory));
addpath(intlabRoot,'-begin');
assert(strcmp(which('startintlab'),fullfile(intlabRoot,'startintlab.m')), ...
    'intlabStartup:Shadow','The selected startintlab is outside the requested root.');
invokeFullStartup();
report = intlab_runtime_check();
end

function invokeFullStartup()
% Isolate caller-workspace variables if a local installation uses a script.
startintlab;
end
