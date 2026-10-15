function ai_middle_pool()
runtime=fileparts(mfilename('fullpath'));addpath(runtime,'-begin');clear chebmag;
if isempty(gcp('nocreate')),ge_verified_pool();end
assert(strcmp(which('chebmag'),fullfile(runtime,'chebmag.m')));
end
