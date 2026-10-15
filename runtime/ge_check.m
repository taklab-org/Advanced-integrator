function report=ge_check()
report.runtime=intlab_runtime_check();
failed=intvalinit('checkrounding',0);
assert(~failed,'GE:Rounding','INTLAB rounding self-test failed.');
report.blas=intlab_blas_check(512);
report.blasVariational=intlab_blas_check(2352);
report.intervalMatrixPassed=report.blas.passed && report.blasVariational.passed;
report.pid=feature('getpid');
report.threads=maxNumCompThreads;
setround(0);
end
