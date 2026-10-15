function report = intlab_blas_check(n)
%INTLAB_BLAS_CHECK Check a nonrepresentable exact matrix product after startup.
% Run in each process at its intended BLAS thread setting. No cache writes.
if nargin==0, n=512; end
validateattributes(n,{'double'},{'scalar','integer','>=',2,'finite'});
patchReport=[];
if strcmp(getenv('TAKLAB_MATLAB_BLAS_POLICY'),'accelerate-rounding-v1')
    patchReport=intlab_accelerate_check();
end
assert(~isempty(which('intval')) && ~isempty(which('setround')), ...
    'intlabStartup:Path','Initialize INTLAB before this matrix probe.');
originalRounding=getround;
restoreRounding=onCleanup(@() setround(originalRounding));
setround(0);
A=zeros(n); B=zeros(n);
A(:,1)=1; B(1,:)=1; A(:,2)=2^-27; B(2,:)=2^-27;
% Each exact entry is 1+2^-54. Its binary64 enclosure is [1,1+eps].
expectedLower=1; expectedUpper=1+eps;
setround(-1); D=A*B;
setround(1); U=A*B;
setround(0); C=intval(A)*B;
lo=inf(C); hi=sup(C);
report=struct('blas',version('-blas'),'matlabVersion',version, ...
    'computer',computer,'threadLimit',maxNumCompThreads,'size',n, ...
    'badDown',nnz(~isfinite(D)|D>expectedLower), ...
    'badUp',nnz(~isfinite(U)|U<expectedUpper), ...
    'badEnclosures',nnz(~isfinite(lo)|~isfinite(hi)|lo>expectedLower|hi<expectedUpper));
report.passed=report.badDown==0 && report.badUp==0 && report.badEnclosures==0;
report.patch=patchReport;
assert(report.passed,'intlabStartup:BlasRounding', ...
    'BLAS matrix rounding failed: %s; n=%d; down=%d up=%d interval=%d.', ...
    report.blas,n,report.badDown,report.badUp,report.badEnclosures);
end
