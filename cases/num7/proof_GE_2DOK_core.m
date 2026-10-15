function proof_GE_2DOK_core(ic_file, N_fixed, phase2_file, phase3_file, result_file, label)
% PROOF_GE_2DOK_CORE  Computer-assisted proof of global existence for the
% 2D Ohta--Kawasaki equation, starting from a given initial condition.
%
%   proof_GE_2DOK_core(ic_file, N_fixed, phase2_file, phase3_file, result_file, label)
%
% Pipeline (see AT-notes/main.tex, Sections 3-4):
%   Phase 1  ODE time integration of the Galerkin system (ode113).
%   Phase 2  Newton-refined Chebyshev approximation + defect / semigroup data.
%   Phase 3  Rigorous variational equation on every time step (parfeval).
%   Phase 4  Forward non-autonomous CAP on [0, tau_F] (bare product + direct
%            recurrence), with K_F chosen automatically as the largest
%            successful truncation.
%   Step 4   Validated intermediate trajectory with nonlinear step tests,
%            retaining initial correlations and finite/tail error bounds.
%   Step 5   Infinite step on [tau_M, inf): radii polynomial at u_*.
%
% If phase2_file AND phase3_file already exist they are loaded; otherwise
% Phase 1-3 are recomputed from ic_file (Phase 3 is the bottleneck, ~days).
% Phase 4, Step 4 and Step 5 are always evaluated in place.
%
% Inputs
%   ic_file      initial-condition .mat (fields: params, N, L, a0_full)
%   N_fixed      Fourier truncation for the fixed-N Phase 2 recomputation
%   phase2_file  Phase 2 data file (load if present, else written)
%   phase3_file  Phase 3 data file (load if present, else written)
%   result_file  output .mat with the verified GE enclosure
%   label        short text label for log output

%% ====================== environment ======================
basedir = fileparts(mfilename('fullpath'));
% Initialized by run_case; all paths and reports are run-local.
ai_add_case_paths(basedir);
fprintf('\n############  GE proof: 2D Ohta--Kawasaki  [%s]  ############\n', label);

%% =============== Phase 1-3: load or compute ===============
if exist(phase2_file,'file') && exist(phase3_file,'file')
  fprintf('Phase 2/3 data present -- loading.\n  %s\n  %s\n', phase2_file, phase3_file);
else
  fprintf('Phase 2/3 data absent -- recomputing from %s.\n', ic_file);
  phaseClock=tic;compute_phase2(ic_file, N_fixed, phase2_file);phase12_seconds=toc(phaseClock);save(fullfile(getenv('GE_VERIFIED_RUN_DIR'),'phase12_timing.mat'),'phase12_seconds');
  phaseClock=tic;compute_phase3(phase2_file, phase3_file);phase3_seconds=toc(phaseClock);save(fullfile(getenv('GE_VERIFIED_RUN_DIR'),'phase3_timing.mat'),'phase3_seconds');
end
[source_data_revalidated,source_certificate]=ge_audit_source_data( ...
  ic_file,N_fixed,phase2_file,phase3_file,basedir);

if strcmp(getenv('AI_SMOKE'),'1')
 assert(source_data_revalidated); save(result_file,'source_certificate','source_data_revalidated'); return;
end
load(phase2_file);                                                       %#ok<LOAD>
P3 = load(phase3_file, 'bUn_vec', 'bUn_endpt_vec', 'r1_vec');
bUn_vec = P3.bUn_vec; bUn_endpt_vec = P3.bUn_endpt_vec; r1_vec = P3.r1_vec;

lambda = params(1); sigma = params(2); epsilon = params(3);
N1 = N(1); N2 = N(2); L1 = L(1); L2 = L(2);
m1_F = m(1); m2_F = m(2);
N_eq = [max(24,N_fixed),max(24,N_fixed)]; % Equilibrium approximation; independent of time/linear blocks.
m_inf = N; mN = prod(m_inf);
assert(all(m_inf<=N_eq),'The linear block must fit the equilibrium coefficient support.');
K_total = num_integration;
fprintf('  m=[%d,%d] N=[%d,%d] K_total=%d tau_total=%.4f epsilon=%g\n', ...
  m1_F, m2_F, N1, N2, K_total, sum(mid(tau_vec)), epsilon);

%% =================== Step 1: equilibrium ==================
fprintf('\n=== Step 1: equilibrium u_* (Newton) ===\n');
N_last = round(sqrt(numel(a_end_final)));
a_end_mat = reshape(a_end_final, N_last, N_last);
a_seed=zeros(N_eq);overlap=min(N_eq,[N_last,N_last]);
a_seed(1:overlap(1),1:overlap(2))=a_end_mat(1:overlap(1),1:overlap(2));
a_ss=a_seed(:);a_ss=a_ss(2:end);
for it = 1:200
  Fv = F_OK_2D(a_ss, params, N_eq, L); nF = norm(Fv, inf);
  if nF < 1e-14, break; end
  DF = DF_ss_OK(a_ss, params, N_eq, L);
  da = DF \ Fv; al = 1;
  for ls = 1:25
    if norm(F_OK_2D(a_ss-al*da,params,N_eq,L), inf) < nF, break; end
    al = al*0.5;
  end
  a_ss = a_ss - al*da;
end
fprintf('  Newton: iter=%d  ||F(u_*)||=%.2e\n', it, norm(F_OK_2D(a_ss,params,N_eq,L),inf));
a_star_full = reshape([lambda; a_ss], N_eq);
ia_star = intval(a_star_full);
[~, ia2_star] = powerconvcos(ia_star, 2);
[~, ia3_full] = powerconvcos(ia_star, 3);
[ikx, iky] = ndgrid(0:N_eq(1)-1, 0:N_eq(2)-1);
ibkL2 = (ikx*intval(L1)).^2 + (iky*intval(L2)).^2;
imu_k  = ibkL2 .* (-intval(epsilon)^2*ibkL2 + 1) - intval(sigma);
% The cubic creates modes outside F_N; include their residual as well.
[nr1, nr2] = size(ia3_full);
[ikx_full, iky_full] = ndgrid(0:nr1-1, 0:nr2-1);
ibkL2_full = (ikx_full*intval(L1)).^2 + (iky_full*intval(L2)).^2;
res_star = -ibkL2_full .* ia3_full;
res_star(1:N_eq(1),1:N_eq(2)) = res_star(1:N_eq(1),1:N_eq(2)) + imu_k.*ia_star;
Y1_inf      = wnorm(res_star, nu);
banorm_star = wnorm(ia_star, nu);
fprintf('  Y1(defect at u_*) = %.4e   ||u_*||_nu = %.4e\n', sup(Y1_inf), sup(banorm_star));

%% ============ Step 2: forward Phase 4, auto K_F ===========
fprintf('\n=== Step 2: forward Phase 4, auto-search K_F ===\n');
delete(gcp('nocreate'))
ge_verified_pool(4, 4)
assert(K_total==128 && sum(tau_vec)==8,'GE:RequestedBoundary','Requested boundary must be step128/time8.'); K_F_list = 128;
scan = struct('K_F',{},'ok',{},'sr',{},'re',{},'verification',{});
fut(1:numel(K_F_list)) = parallel.FevalFuture;
for k = 1:numel(K_F_list)
  fut(k) = parfeval(@phase4_worker, 7, K_F_list(k), bUn_vec, bUn_endpt_vec, r1_vec, ...
    int_semi_store, int_semi_q_store, psi_vec_store, c_out_store, g22_store, ...
    VinfJ_vec, VinfqJ_vec, F1ba_fin_vec, F1ba_inf_vec, F2ba_fin_vec, F2ba_inf_vec, ...
    banorm_vec, tau_vec, mu_m, m, N, nu);
end
for c = 1:numel(K_F_list)
  [kl, ok, sr, ~, ~, ~, re, verification] = fetchNext(fut);
  scan(end+1) = struct('K_F',K_F_list(kl),'ok',ok,'sr',sr,'re',re,'verification',verification);   %#ok<AGROW>
  save([result_file '.scan.mat'],'scan');
  fprintf('  K_F=%d  ok=%d  seg_r=%.3e  r_end=%.3e\n', K_F_list(kl), ok, sr, re);
end
okv = [scan.K_F]; okm = [scan.ok]==1;
if ~any(okm), delete(gcp('nocreate')); error('no K_F succeeded in coarse scan'); end
K_F_lo = max(okv(okm));
failhi = okv(~okm & okv>K_F_lo);
if ~isempty(failhi)
  K_ref = (min(failhi)-1):-1:(K_F_lo+1);
  fut2(1:numel(K_ref)) = parallel.FevalFuture;
  for k = 1:numel(K_ref)
    fut2(k) = parfeval(@phase4_worker, 7, K_ref(k), bUn_vec, bUn_endpt_vec, r1_vec, ...
      int_semi_store, int_semi_q_store, psi_vec_store, c_out_store, g22_store, ...
      VinfJ_vec, VinfqJ_vec, F1ba_fin_vec, F1ba_inf_vec, F2ba_fin_vec, F2ba_inf_vec, ...
      banorm_vec, tau_vec, mu_m, m, N, nu);
  end
  for c = 1:numel(K_ref)
    [kl, ok, sr, ~, ~, ~, re, verification] = fetchNext(fut2);
    scan(end+1) = struct('K_F',K_ref(kl),'ok',ok,'sr',sr,'re',re,'verification',verification);     %#ok<AGROW>
    save([result_file '.scan.mat'],'scan');
    fprintf('  K_F=%d  ok=%d  seg_r=%.3e  r_end=%.3e  (refine)\n', K_ref(kl), ok, sr, re);
  end
end
delete(gcp('nocreate'))
K_F = max([scan([scan.ok]==1).K_F]);
ixF = find([scan.K_F]==K_F & [scan.ok]==1, 1);
re_F_poly=scan(ixF).re;sr_F=scan(ixF).sr;

forward_verification=scan(ixF).verification;
tau_F_bound=sum(intval(tau_vec(1:K_F)));tau_F=sup(tau_F_bound);
fprintf('  selected K_F=%d  tau_F=%.4f  r_end=%.4e\n', K_F, tau_F, re_F_poly);
assert(forward_verification.radii_verified && forward_verification.composition_verified);phase4_source=ge_source_manifest(basedir);save(fullfile(getenv('GE_VERIFIED_RUN_DIR'),'phase4_selected.mat'),'K_F','tau_F','re_F_poly','forward_verification','scan','phase4_source','-v7.3');
fid=fopen(fullfile(getenv('GE_VERIFIED_RUN_DIR'),'PHASE4_PASSED'),'w');fprintf(fid,'Actual forward radius and composition selected K=%d.\n',K_F);fclose(fid);

% matching error at tau_F
N_F1 = round(sqrt(numel(a0_init_vec{K_F+1})));
a_F_mat = reshape(a0_init_vec{K_F+1}, N_F1, N_F1);
% The CAP radius is about the exact polynomial trace; the next candidate
% starts from a nominal floating endpoint. Enclose this additional gap.
last_ba=intval(ba_vec{K_F});trace_N=[size(last_ba,2),size(last_ba,3)];
trace_end=reshape(last_ba(1,:,:)+2*sum(last_ba(2:end,:,:),1),trace_N);
common=max(trace_N,[N_F1,N_F1]);gap_F=intval(zeros(common));
gap_F(1:trace_N(1),1:trace_N(2))=trace_end;
gap_F(1:N_F1,1:N_F1)=gap_F(1:N_F1,1:N_F1)-intval(a_F_mat);
endpoint_gap_F=wnorm(gap_F,intval(nu));re_F=sup(intval(re_F_poly)+endpoint_gap_F);
N_match=max(N_eq,[N_F1,N_F1]);
ia_star_pad=intval(zeros(N_match));ia_star_pad(1:N_eq(1),1:N_eq(2))=ia_star;
ia_F_pad=intval(zeros(N_match));ia_F_pad(1:N_F1,1:N_F1)=intval(a_F_mat);
Y2_F = wnorm(ia_star_pad-ia_F_pad,intval(nu));
rho_0 = intval(re_F) + Y2_F;
fprintf('  Y2(||u_* - a(tau_F)||) = %.4e   rho_0 = %.4e\n', sup(Y2_F), sup(rho_0));

%% ======== Step 3: parameters of the validated intermediate phase =======
m_a=[12,12];h_auto=0.1;K_auto=6000;
ba2_sup_star=abs(ia2_star);
% The reference trajectory keeps the forward endpoint's Fourier support.
% Its omitted modes are covered by the FULL-PDE defect and tail error bounds.
assert(all(m_a<=[N_F1,N_F1]));

%% ===== Precompute terminal bounds before choosing a bridge endpoint =====
fprintf('\n=== Terminal bounds (m_inf=[%d,%d]) ===\n', m_inf(1), m_inf(2));
[kxN,kyN] = ndgrid(0:m_inf(1)-1,0:m_inf(2)-1); kxvN=kxN(:); kyvN=kyN(:);
alpN=4*ones(m_inf(1),m_inf(2)); alpN(1,:)=2; alpN(:,1)=2; alpN(1,1)=1;
wN=alpN(:).*nu.^kxvN.*nu.^kyvN; D_wN=diag(intval(wN)); D_wN_inv=diag(1./intval(wN));
ibkL2_N=(kxN*intval(L1)).^2+(kyN*intval(L2)).^2;
imu_N=ibkL2_N.*(-intval(epsilon)^2*ibkL2_N+1)-intval(sigma);
DN_N=intval(zeros(mN));
for j=1:mN, l1=kxvN(j); l2=kyvN(j);
  for i=1:mN, k1=kxvN(i); k2=kyvN(i);
    DN_N(i,j) = DN_entry(ia2_star, k1,k2,l1,l2);
  end
end
AN=diag(imu_N(:))-ibkL2_N(:).*DN_N;
[VN,DNe]=eig(mid(AN)); ilamN=intval(diag(DNe)); iVN=intval(VN); iVNi=inv(iVN);
if ~all(sup(real(ilamN)) < 0)
  error('Step 5: finite-block eigenvalues are not strictly stable');
end
[VinfJ_N,VinfqJ_N]=getting_Vinf(Inf,max(2*m_inf,m_inf+1),m_inf,params,L);
inv_lamN=-1./real(ilamN);
V_tilN = abs(iVN)*diag(inv_lamN);
M_N    = V_tilN*(abs(iVNi)*abs(AN-iVN*diag(ilamN)*iVNi));   % M1 defect matrix (>=0)
B_neuN=intval(eye(mN)) - M_N;
Mn_N = sup(norm(M_N,inf));
if Mn_N>=1
  error('Step 5: B_neuN not an M-matrix / (I-M)^{-1} not nonnegative: ||M||_inf=%.4f >= 1', Mn_N);
end
int_semiN=B_neuN\(V_tilN*abs(iVNi));
int_semi_wN =D_wN*int_semiN*D_wN_inv;
int_semi_wqN=D_wN*(int_semiN.*ibkL2_N(:).')*D_wN_inv;
% ---- M1 rigorous Phi_sup on [0,inf): interval-subinterval + infinite tail block ----
relam_N = real(ilamN);
% geometric grid: first node t1 s.t. max|Re lam|*t1 = 0.05 (uniform <=5% relative
% width of e^{Re lam*T_i} on every subinterval), ratio 1.1 up to T_max=100
t1_N = 0.05/max(sup(abs(relam_N))); ratio_N = 1.1; Tmax_N = 100;   % double grid nodes
nsub_N = ceil(log(Tmax_N/t1_N)/log(ratio_N));
tsN = [0, t1_N*ratio_N.^(0:nsub_N)];
Sbar_N = zeros(mN);
absV_N = abs(iVN); absVi_N = abs(iVNi);
for ig=1:numel(tsN)-1
  % 2nd-order Taylor form: sup_{t in [ta,tb]} |C(t)|
  %   <= |C(tm)| + (D/2)|C'(tm)| + (D^2/8) sup_{[ta,tb]}|C''|,   D = tb-ta,
  % with C'(tm)=V*diag(lam.*e^{lam tm})*V^{-1} SIGNED (keeps cancellation) and
  % |C''(t)| <= |V|*diag(lam.^2.*e^{lam t})*|V^{-1}| majorized at t=ta
  % (lam^2 e^{lam t} decreasing for lam<0). Only the 2nd-order term pays |V||V^{-1}|.
  ta = intval(tsN(ig)); tb = intval(tsN(ig+1)); tm = (ta+tb)/2; D = tb-ta;
  Cmid  = abs(iVN*diag(exp(relam_N*tm))*iVNi);
  Cpmid = abs(iVN*diag(relam_N.*exp(relam_N*tm))*iVNi);
  d2bnd = absV_N*diag(relam_N.^2.*exp(relam_N*ta))*absVi_N;
  Sbar_N = max(Sbar_N, sup(Cmid + (D/2).*Cpmid + (D^2/8).*d2bnd));
end
% tail t in [tsN(end),inf): e^{Re lam t} in [0, e^{Re lam*tsN(end)}] (Re lam<0)
etN = infsup(zeros(mN,1), sup(exp(relam_N*intval(tsN(end)))));
Sbar_N = max(Sbar_N, sup(abs(iVN*diag(etN)*iVNi)));
Sbar_N = intval(Sbar_N);
Phi_sup_N = B_neuN\Sbar_N;                       % S <= B_neuN^{-1} Sbar
Phi_sup_wN=D_wN*Phi_sup_N*D_wN_inv;
fprintf('  [M1 block N] ||M||_inf=%.4e  ||Phi_sup_wN||_1=%.4e\n', ...
  Mn_N, sup(max(sum(abs(Phi_sup_wN),1))));
[psi_N,cout_N]=getting_coupling_vectors(ia2_star, m_inf, N_eq, nu, 'signed');
cin_wN =D_wN*(abs(int_semiN.*ibkL2_N(:).')*(3*psi_N));
cout_wN=(D_wN_inv*(VinfqJ_N*3*cout_N))';
schur_N=intval(1)-VinfqJ_N*3*wnorm(ba2_sup_star,nu)-cout_wN*cin_wN;
if ~(inf(schur_N) > 0)
  error('Step 5: infinite-step Schur complement is not strictly positive');
end
sinv_N=1/schur_N;
BN=intval(zeros(mN+1));
BN(1:mN,1:mN)=intval(eye(mN))+cin_wN*(sinv_N*cout_wN);
BN(1:mN,mN+1)=cin_wN*sinv_N; BN(mN+1,1:mN)=sinv_N*cout_wN; BN(mN+1,mN+1)=sinv_N;
bUN =intval(zeros(mN+1)); bUN(1:mN,1:mN)=Phi_sup_wN;    bUN(mN+1,mN+1)=intval(1);
bVqN=intval(zeros(mN+1)); bVqN(1:mN,1:mN)=int_semi_wqN; bVqN(mN+1,mN+1)=VinfqJ_N;
bV0N=intval(zeros(mN+1)); bV0N(1:mN,1:mN)=int_semi_wN;  bV0N(mN+1,mN+1)=VinfJ_N;
W_inf =sup(max(sum(abs(BN*bUN ),1)));
Vq_inf=sup(max(sum(abs(BN*bVqN),1)));
V0_inf=sup(max(sum(abs(BN*bV0N),1)));
fprintf('  W_inf=%.4e  V0_inf=%.4e  Vq_inf=%.4e  Schur=%.4e\n', ...
  W_inf, V0_inf, Vq_inf, inf(schur_N));

forcing_inf=getting_constant_forcing(res_star,m_inf,Inf,params,L,nu,BN,int_semi_wN);
Y_f_inf=sum(forcing_inf);
fprintf('  mode-resolved residual Volterra bound=%.4e (raw residual norm %.4e)\n', ...
  sup(Y_f_inf),sup(Y1_inf));

%% ============= Step 4: sequential nonlinear intermediate CAP ==========
terminal=struct('W',intval(W_inf),'forcing',Y_f_inf, ...
  'Vq',intval(Vq_inf),'anorm',banorm_star);
save([result_file '.forward.mat'],'label','K_F','tau_F','re_F','sr_F','Y2_F','rho_0', ...
  're_F_poly','endpoint_gap_F','tau_F_bound','forward_verification','source_certificate', ...
  'source_data_revalidated','a_F_mat','a_star_full','params','L','nu','m_a','h_auto','K_auto', ...
  'terminal','N_eq','m_inf','W_inf','V0_inf','Vq_inf','Y_f_inf','scan','-v7.3');
bridge_options=struct('progress_file',[result_file '.bridge_progress.mat']);
ai_middle_pool();
bridge=verify_intermediate_chain(a_F_mat,intval(re_F),a_star_full,params,L,nu, ...
  m_a,h_auto,K_auto,terminal,bridge_options);
K_bridge=bridge.steps;rho_bridge=sup(bridge.rho);bridge_test=bridge.terminal_test;
bridge_history=bridge.history;

%% ============= Step 5: terminal local test and proof status =============
Y_infty=bridge_test.Y;Pr=bridge_test.Pr;Zv=bridge_test.Zv;r_up=bridge_test.radius;
tau_total=sup(tau_F_bound+bridge.duration);
validation_status=struct('forward_radii_passed',forward_verification.radii_verified, ...
  'forward_composition_verified',forward_verification.composition_verified, ...
  'bridge_nonlinear_verified',bridge.verified,'terminal_radii_passed',bridge_test.passed, ...
  'source_data_revalidated',source_data_revalidated);
% Legacy caches remain uncertified. Fresh files require matching generation
% records, initial data, sources, contraction evidence, and runtime checks.
success=validation_status.forward_radii_passed ...
  && validation_status.forward_composition_verified ...
  && validation_status.bridge_nonlinear_verified && validation_status.terminal_radii_passed ...
  && validation_status.source_data_revalidated;
save(result_file, 'label','K_F','tau_F','re_F','sr_F','Y2_F','rho_0', ...
  'm_a','K_bridge','rho_bridge','tau_total','W_inf','V0_inf','Vq_inf', ...
  'Pr','Zv','r_up','a_star_full','N_eq','m_inf','Y_f_inf','Y_infty', ...
  'bridge_test','bridge_history','bridge','forward_verification', ...
  're_F_poly','endpoint_gap_F','tau_F_bound','validation_status','source_certificate','success','scan','-v7.3');

fprintf('\n=================================================================\n');
if success
  fprintf('  ===  GLOBAL EXISTENCE PROOF SUCCESS  [%s]  ===\n', label);
  fprintf('  Forward [0, %.4f]      K_F=%d, r_end=%.4e\n', tau_F, K_F, re_F);
  fprintf('  Intermediate [%.4f, %.4f]  %d steps, m_a=[%d,%d]\n', ...
    tau_F, tau_total, K_bridge, m_a(1), m_a(2));
  fprintf('  Step 5  [%.4f, inf)    r_inf=%.4e, Z=%.4e < 1\n', tau_total, r_up, sup(Zv));
  fprintf('  ||u(t) - u_*||_nu <= %.4e   for all t >= %.4f\n', r_up, tau_total);
else
  fprintf('  ===  GLOBAL PROOF NOT VERIFIED  [%s]  ===\n', label);
  fprintf('  Terminal local test: P=%.4e, Z=%.4e, radius=%.4e\n',sup(Pr),sup(Zv),r_up);
  fprintf('  Intermediate validation: %d, terminal test: %d\n',bridge.verified,bridge_test.passed);
  fprintf('  Source-data revalidation: %d\n',source_data_revalidated);
  if ~bridge.verified,fprintf('  %s\n',bridge.reason);end
end
fprintf('  result saved to %s\n', result_file);
fprintf('=================================================================\n');
end

%% ===================================================================
%%                         local functions
%% ===================================================================
function compute_phase2(ic_file, N_fixed, phase2_file)
% Phase 1 (ODE) + Phase 2 (Newton refinement, defect, semigroup data),
% fixed Fourier truncation N_fixed, fixed time step tau1 = 2^-4, target tau = 13.
  S = load(ic_file);
  params = S.params; L = S.L; nu = 1.0;
  m = [7,7]; m1=m(1); m2=m(2);
  tau = 8; if strcmp(getenv('AI_SMOKE'),'1'),tau=2^-4;end; tau1 = 2^-4;
  rel_tol = 1e-14; abs_floor = 1e-13;
  N_init = S.N(1);
  fprintf('  Phase 1-2: N=%d, m=[7,7], tau=%.4f, tau1=%.4f\n', N_fixed, tau, tau1);

  [kxw,kyw] = ndgrid(0:m1-1,0:m2-1);
  alp = 4*ones(m1,m2); alp(1,:)=2; alp(:,1)=2; alp(1,1)=1;
  w_fin = alp(:).*nu.^kxw(:).*nu.^kyw(:);
  D_w = diag(intval(w_fin)); D_w_inv = diag(1./intval(w_fin));

  Kmax = 250;
  ba_vec=cell(Kmax,1); a0_init_vec=cell(Kmax+1,1);
  banorm_vec=intval(zeros(Kmax,1)); tau_vec=zeros(Kmax,1);
  F1ba_fin_vec=intval(zeros(Kmax,1)); F1ba_inf_vec=intval(zeros(Kmax,1));
  F2ba_fin_vec=intval(zeros(Kmax,1)); F2ba_inf_vec=intval(zeros(Kmax,1));
  VinfJ_vec=intval(zeros(Kmax,1));    VinfqJ_vec=intval(zeros(Kmax,1));
  int_semi_store=cell(Kmax,1); int_semi_q_store=cell(Kmax,1);
  psi_vec_store=cell(Kmax,1);  c_out_store=cell(Kmax,1);
  g22_store=intval(zeros(Kmax,1));
  N_per_seg=zeros(Kmax,1); n_cheb_per_seg=zeros(Kmax,1);

  a0_pad = zeros(N_fixed,N_fixed);
  a0_init = reshape(S.a0_full(:), N_init, N_init);
  nmin = min(N_init, N_fixed);
  a0_pad(1:nmin,1:nmin) = a0_init(1:nmin,1:nmin);
  a0_init_vec{1} = a0_pad(:);

  previous_endpoint_error=intval(zeros(N_fixed));
  source_at_start=ge_source_manifest(fileparts(mfilename('fullpath')));
  checkpoint_dir=[phase2_file '.steps'];if ~exist(checkpoint_dir,'dir'),mkdir(checkpoint_dir);end
  schur_margins=zeros(Kmax,1);
  delete(gcp('nocreate'));
  ge_verified_pool(4,4);
  futures=parallel.FevalFuture.empty(0,1);
  mu_m_max=-inf;ts=1;tstart=0;t0=tic;
  while tstart<tau && ts<=Kmax
    tau1u=min(tau1,tau-tstart);
    candidate=build_segment_candidate([N_fixed,N_fixed],[tstart,tstart+tau1u], ...
      a0_init_vec{ts},params,L,rel_tol,abs_floor);
    a0_init_vec{ts+1}=candidate.a_end;tau_vec(ts)=tau1u;
    futures(ts)=parfeval(@verify_segment_candidate,2,candidate,params,L,m,nu,D_w,D_w_inv,previous_endpoint_error);
    previous_endpoint_error=candidate.end_error;
    if mod(ts,10)==0,fprintf('    Candidate %d ready, tau=%.4f; verification dispatched\n',ts,tstart+tau1u);end
    tstart=tstart+tau1u;ts=ts+1;
  end
  K=ts-1;
  for completed=1:K
    [j,seg,F1]=fetchNext(futures);ge_validate_segment(seg);schur_margins(j)=seg.schur_margin;
    tau1u=tau_vec(j);segment_start=sum(tau_vec(1:j-1));
    save(fullfile(checkpoint_dir,sprintf('step_%04d.mat',j)),'seg','j','segment_start','tau1u','params','L','nu','m','-v7.3');
    ba_vec{j}=seg.ba;banorm_vec(j)=seg.banorm;
    F1ba_fin_vec(j)=seg.F1_fin;F1ba_inf_vec(j)=seg.F1_inf;
    F2ba_fin_vec(j)=seg.F2_fin;F2ba_inf_vec(j)=seg.F2_inf;
    VinfJ_vec(j)=seg.VinfJ;VinfqJ_vec(j)=seg.VinfqJ;
    int_semi_store{j}=seg.int_semi;int_semi_q_store{j}=seg.int_semi_q;
    psi_vec_store{j}=seg.psi_vec;c_out_store{j}=seg.c_out_pm;g22_store(j)=seg.g22;
    N_per_seg(j)=N_fixed;n_cheb_per_seg(j)=seg.n_cheb;
    mu_m_max=max(mu_m_max,sup(seg.mu_m_seg));
    fprintf('    Phase 2 verified %d/%d (segment %d): n_cheb=%d F1=%.3e elapsed=%.1fs\n',completed,K,j,seg.n_cheb,F1,toc(t0));
  end
  delete(gcp('nocreate'));
  num_integration = ts-1;
  idx = 1:num_integration;
  ba_vec=ba_vec(idx); banorm_vec=banorm_vec(idx); tau_vec=tau_vec(idx);
  F1ba_fin_vec=F1ba_fin_vec(idx); F1ba_inf_vec=F1ba_inf_vec(idx);
  F2ba_fin_vec=F2ba_fin_vec(idx); F2ba_inf_vec=F2ba_inf_vec(idx);
  VinfJ_vec=VinfJ_vec(idx); VinfqJ_vec=VinfqJ_vec(idx);
  int_semi_store=int_semi_store(idx); int_semi_q_store=int_semi_q_store(idx);
  psi_vec_store=psi_vec_store(idx); c_out_store=c_out_store(idx);
  g22_store=g22_store(idx); N_per_seg=N_per_seg(idx); n_cheb_per_seg=n_cheb_per_seg(idx);
  a0_init_vec=a0_init_vec(1:num_integration+1);
  t_cumul=cumsum(tau_vec); mu_m=intval(mu_m_max);
  a_end_final=a0_init_vec{num_integration+1}; N=[N_fixed,N_fixed];
  assert(tstart==tau && num_integration>0,'GE:Phase2Incomplete','The requested time interval was not completed.');
  assert(isequal(source_at_start,ge_source_manifest(fileparts(mfilename('fullpath')))),'GE:SourceChanged','Source changed during Phase 2.');
  generation_record=struct('schema',1,'phase',2,'complete',true,'source',source_at_start, ...
    'ic_sha256',ge_file_sha256(ic_file),'N',N,'m',m,'steps',num_integration, ...
    'schur_margins',schur_margins(idx),'target_time',tau);
  fprintf('  Phase 2 done in %.1f s, K=%d\n', toc(t0), num_integration);
  save(phase2_file,'num_integration','ba_vec','banorm_vec','tau_vec', ...
    'F2ba_fin_vec','F2ba_inf_vec','F1ba_fin_vec','F1ba_inf_vec', ...
    'VinfJ_vec','VinfqJ_vec','int_semi_store','int_semi_q_store', ...
    'psi_vec_store','c_out_store','g22_store','m','N','nu','params','L', ...
    'mu_m','t_cumul','tau1','a0_init_vec','a_end_final','N_per_seg', ...
    'n_cheb_per_seg','N_fixed','generation_record','-v7.3');
end

function compute_phase3(phase2_file, phase3_file)
% Phase 3: rigorous variational equation on every time step via parfeval.
  P2 = load(phase2_file,'ba_vec','tau_vec','m','params','L','num_integration');
  phase2_sha256=ge_file_sha256(phase2_file);
  source_at_start=ge_source_manifest(fileparts(mfilename('fullpath')));
  checkpoint_dir=[phase3_file '.steps'];if ~exist(checkpoint_dir,'dir'),mkdir(checkpoint_dir);end
  m1m2 = prod(P2.m); K = P2.num_integration;
  fprintf('  Phase 3: %d variational steps (parfeval 4x4)\n', K);
  delete(gcp('nocreate'))
  ge_verified_pool(4,4)
  bUn_vec=cell(K,1); bUn_endpt_vec=cell(K,1); r1_vec=intval(zeros(K,m1m2));
  verification_vec=cell(K,1);
  fut(1:K) = parallel.FevalFuture;
  for ts = 1:K
    fut(ts) = parfeval(@solve_variational_step, 4, ...
      P2.ba_vec{ts}, P2.tau_vec(ts), P2.m, P2.params, P2.L);
  end
  t0 = tic;
  for c = 1:K
    [ts, bUn, bUn_e, r1, verification] = fetchNext(fut);
    assert(verification.passed && all(inf(verification.margin)>0,'all'));
    verification_vec{ts}=verification;
    save(fullfile(checkpoint_dir,sprintf('step_%04d.mat',ts)),'ts','bUn','bUn_e','r1','verification','-v7.3');
    bUn_vec{ts}=bUn; bUn_endpt_vec{ts}=bUn_e; r1_vec(ts,:)=r1;
    if mod(c,20)==0
      fprintf('    [%6.0fs] %d/%d done\n', toc(t0), c, K);
    end
  end
  delete(gcp('nocreate'))
  num_integration = K;                                                   %#ok<NASGU>
  assert(isequal(source_at_start,ge_source_manifest(fileparts(mfilename('fullpath')))));
  assert(strcmp(phase2_sha256,ge_file_sha256(phase2_file)));
  generation_record=struct('schema',1,'phase',3,'complete',true,'source',source_at_start, ...
    'phase2_sha256',phase2_sha256,'steps',K,'verification',{verification_vec});
  save(phase3_file,'bUn_vec','bUn_endpt_vec','r1_vec','num_integration','generation_record','-v7.3');
  fprintf('  Phase 3 done in %.1f s\n', toc(t0));
end

function candidate = build_segment_candidate(N,tspan,a0,params,L,rel_tol,abs_floor)
% One Phase 1+2 segment: ODE solve, Chebyshev projection, Newton refinement,
% defect (F1, F2) and integrated-semigroup / coupling data.
  N1=N(1); N2=N(2); h=tspan(2)-tspan(1);
  lambda=params(1);
  opts = odeset('abstol',1e-16,'reltol',2.22045e-14);
  Y = ode113(@(t,y) F_OK_2D(y,params,N,mid(L)), tspan, a0(2:end), opts);
  a = pdecheb2_relative_v3(Y, tspan, rel_tol, abs_floor);
  n = size(a,1);
  a_full = zeros(n, N1*N2); a_full(1,1)=lambda; a_full(:,2:end)=a;
  [a_full, a_end] = refine_newton(a_full.', a0, n, N, L, params, h);
  a_full = a_full.'; a_iv = intval(a_full);
  ta=[a_iv(1,:);a_iv(2:end,:)/2];ba=reshape(ta,n,N1,N2);
  candidate=struct('N',N,'h',h,'a0',a0,'a_iv',a_iv,'a_end',a_end,'n',n, ...
    'end_error',reshape(ba(1,:,:)+2*sum(ba(2:end,:,:),1),N)-reshape(intval(a_end),N));
end

function [seg,F1]=verify_segment_candidate(candidate,params,L,m,nu,D_w,D_w_inv,previous_endpoint_error)
  N=candidate.N;N1=N(1);N2=N(2);h=candidate.h;a0=candidate.a0;
  a_iv=candidate.a_iv;a_end=candidate.a_end;n=candidate.n;m1=m(1);m2=m(2);
  du = chebdiff(a_iv)/(h/2);
  du = [du(1,:); du(2:end,:)/2];
  du = reshape(du, n-1, N1, N2);
  du_ext = intval(zeros(3*n-2, 3*N1-2, 3*N2-2));
  du_ext(1:n-1,1:N1,1:N2) = du;
  ta = [a_iv(1,:); a_iv(2:end,:)/2];
  defect = du_ext - F_OK_2D_ext(ta, intval(params), n, N, intval(L));
  F1_fin = wnorm(defect(:,1:m1,1:m2), nu);
  defect(:,1:m1,1:m2) = 0;
  F1_inf = wnorm(defect, nu);
  F1 = sup(F1_fin + F1_inf);
  ba = reshape(ta, n, N1, N2);
  veps = compute_eps2(intval(ba), intval(a0))-previous_endpoint_error;
  F2_fin = wnorm(veps(1:m1,1:m2), nu);
  veps(1:m1,1:m2) = 0;
  F2_inf = wnorm(veps, nu);
  banorm = wnorm(reshape(chebmag_oneside(ba),N), nu);
  [int_semi,int_semi_q,ba2_sup] = getting_int_semi_matrix(ba, h, m, params, L, nu);
  [psi_vec,c_out_pm] = getting_coupling_vectors(ba2_sup, m, N, nu);
  [VinfJ,VinfqJ,mu_m_seg] = getting_Vinf(h, N, m, params, L);
  g22 = VinfqJ*3*wnorm(ba2_sup,nu);
  c_in_w  = D_w*(abs(int_semi_q)*(3*psi_vec));
  c_out_w = (D_w_inv*(VinfqJ*3*c_out_pm))';
  schur = intval(1)-g22-c_out_w*c_in_w;
  if ~(inf(schur) > 0)
    error('Phase 2: Schur complement is not strictly positive');
  end
  seg = struct('ba',ba,'a_end',a_end,'banorm',banorm,'F1_fin',F1_fin, ...
    'F1_inf',F1_inf,'F2_fin',F2_fin,'F2_inf',F2_inf,'int_semi',int_semi, ...
    'int_semi_q',int_semi_q,'psi_vec',psi_vec,'c_out_pm',c_out_pm, ...
    'VinfJ',VinfJ,'VinfqJ',VinfqJ,'g22',g22,'mu_m_seg',mu_m_seg,'n_cheb',n, ...
    'end_error',reshape(ba(1,:,:)+2*sum(ba(2:end,:,:),1),N)-reshape(intval(a_end),N), ...
    'schur_margin',inf(schur));
end

function [ok,sr,sZ,Ym,U11,re,verification] = phase4_worker(K_end, bUn_vec, bUn_endpt_vec, r1_vec, ...
    int_semi_store, int_semi_q_store, psi_vec_store, c_out_store, g22_store, ...
    VinfJ_vec, VinfqJ_vec, F1ba_fin_vec, F1ba_inf_vec, F2ba_fin_vec, F2ba_inf_vec, ...
    banorm_vec, tau_vec, mu_m, m, N, nu)
% parfeval wrapper: forward Phase 4 on [0, tau_{K_end}].
  try
    [ok,sr,sZ,Ym,U11,re,verification] = try_phase4(1, K_end, bUn_vec, bUn_endpt_vec, r1_vec, ...
      int_semi_store, int_semi_q_store, psi_vec_store, c_out_store, g22_store, ...
      VinfJ_vec, VinfqJ_vec, F1ba_fin_vec, F1ba_inf_vec, F2ba_fin_vec, F2ba_inf_vec, ...
      banorm_vec, tau_vec, mu_m, m, N, nu, 0);
  catch err
    verification=struct('composition_verified',false,'radii_verified',false,'reason',err.message);
    fprintf(2,'Phase 4 K=%d exception [%s]: %s\n',K_end,err.identifier,err.message);
    ok=false; sr=NaN; sZ=NaN; Ym=NaN; U11=NaN; re=NaN;
  end
end

function v = DN_entry(ia2, k1,k2,l1,l2)
% (i,j) entry of the cosine-convolution Jacobian 3*conv(u_*^2, .) block.
% ia2 contains the full square, including the k+l convolution coefficients.
  if l1==0 && l2==0
    v = 3*ia2(k1+1,k2+1);
  elseif l2==0
    v = 3*(ia2(abs(k1-l1)+1,k2+1) + ia2(k1+l1+1,k2+1));
  elseif l1==0
    v = 3*(ia2(k1+1,abs(k2-l2)+1) + ia2(k1+1,k2+l2+1));
  else
    i1m=abs(k1-l1)+1; i1p=k1+l1+1;
    i2m=abs(k2-l2)+1; i2p=k2+l2+1;
    v = 3*(ia2(i1m,i2m)+ia2(i1p,i2m)+ia2(i1m,i2p)+ia2(i1p,i2p));
  end
end

function DF = DF_ss_OK(a,params,N,L)
% Keep the equilibrium interface; use the same Jacobian as candidate stepping.
DF=DF_OK_2D(a,params,N,L);
end
