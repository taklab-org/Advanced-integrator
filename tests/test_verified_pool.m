function test_verified_pool(out,clientReport)
% Requires a fresh client and the project startup guard selected in MATLABPATH.
if ~isfolder(out),mkdir(out);end
save(fullfile(out,'client_environment.mat'),'clientReport');
setenv('GE_VERIFIED_RUN_DIR',out);
setenv('AI_OPTIONS',jsonencode(struct('Workers',2,'Threads',2,'ClientThreads',2)));
reports=ge_verified_pool();
assert(numel(reports)==2);
for k=1:2
 assert(reports{k}.startup.cacheUnchanged);
 assert(reports{k}.check.intervalMatrixPassed);
end
p=gcp('nocreate');delete(p);
fprintf('VERIFIED_PROCESS_POOL_REGRESSION_PASSED\n');
end
