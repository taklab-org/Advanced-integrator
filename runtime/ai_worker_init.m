function report=ai_worker_init(clientReport,threads)
maxNumCompThreads(threads);
global INTLAB_CONST
assert(isempty(INTLAB_CONST),'A worker startup file initialized INTLAB unexpectedly.');
if isfile(fullfile(clientReport.intlabRoot,'startintlab_minimum.m'))
 startup=intlab_worker_start(clientReport);startup.mode='minimum_prepared_cache';
else
 % Standard INTLAB installations need no nonstandard minimum initializer.
 % ge_verified_pool serializes this full startup across the workers.
 startup=intlab_client_start_worker(clientReport.intlabRoot);
end
check=ge_check();report=struct('startup',startup,'check',check);
end
function r=intlab_client_start_worker(root)
here=pwd;c=onCleanup(@()cd(here));addpath(root,'-begin');invoke();
r=intlab_runtime_check();before=hash(r.cachePath);check=intlab_runtime_check();
r.cacheUnchanged=strcmp(before,hash(check.cachePath));
assert(r.cacheUnchanged);r.mode='serialized_full_startup';
end
function invoke()
startintlab;
end
function h=hash(f)
fid=fopen(f,'rb');assert(fid>0);c=onCleanup(@()fclose(fid));b=fread(fid,Inf,'*uint8');
md=java.security.MessageDigest.getInstance('SHA-256');md.update(typecast(b,'int8'));
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
