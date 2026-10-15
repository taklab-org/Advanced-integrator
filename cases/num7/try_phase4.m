function [ok, seg_r, seg_Z, Y_max, U11_K1, r_end, verification] = try_phase4( ...
    K_start, K_end, ...
    glob_bUn_vec, glob_bUn_endpt_vec, glob_r1_vec, ...
    int_semi_store, int_semi_q_store, psi_vec_store, c_out_store, g22_store, ...
    VinfJ_vec, VinfqJ_vec, ...
    F1ba_fin_vec, F1ba_inf_vec, F2ba_fin_vec, F2ba_inf_vec, ...
    banorm_vec, tau_vec, mu_m, m, N, nu, r_prev)
% Phase 4 for 2D-OK: (m1m2+1) directional resolvent with dual semigroup.
%
% Uses TWO vloc_res matrices per step:
%   vloc_res   - int_semi_mat (for Y-bound, defect already has bkL2)
%   vloc_res_q - int_semi_mat_q (for Z-bound, adds bkL2 for Lipschitz)
%
% All block norms and accumulations are in intval for rigorous bounds.
%
% Returns r_end = r_upper(end) (endpoint radius at step K_end) for
% unified infinite-step verification. 0 if not converged.

if nargin < 23, r_prev = 0; end
verification=struct('composition_verified',false,'radii_verified',false);

m1 = m(1); m2 = m(2); m1m2 = m1*m2;
K = K_end - K_start + 1;

% 2D weight matrices
[kx_w, ky_w] = ndgrid(0:m1-1, 0:m2-1);
alp_2d = 4*ones(m1, m2); alp_2d(1,:) = 2; alp_2d(:,1) = 2; alp_2d(1,1) = 1;
w_fin = alp_2d(:) .* intval(nu).^kx_w(:) .* intval(nu).^ky_w(:);
D_w = diag(intval(w_fin));
D_w_inv = diag(1 ./ intval(w_fin));

%% Step 1: Build local directional resolvents
sup_res = cell(K,1);
vloc_res = cell(K,1); vloc_res_q = cell(K,1);
Phi_endpt_w_store = cell(K,1);
sigma_store = intval(zeros(K,1));
remainder_blocks=intval(zeros(2,2,K));
schur_ok = true;

for ii = 1:K
  if strcmp(getenv('GE_DIAGNOSTICS'),'1') && mod(ii,20)==0,fprintf('phase4 local %d/%d\n',ii,K);end
  ts = K_start + ii - 1;
  h_ts = tau_vec(ts);
  sigma_step = exp(intval(mu_m) * intval(h_ts));

  % Weighted supremum and endpoint of fundamental matrix
  bUn = glob_bUn_vec{ts};
  ps = intval(zeros(m1m2, m1m2));
  for jj = 1:m1m2
    ps(:,jj) = chebmag_oneside(bUn(:,:,jj));
  end
  Phi_sup_w = D_w * (ps + glob_r1_vec(ts,:)) * D_w_inv;
  Phi_endpt_w = D_w * intval(glob_bUn_endpt_vec{ts}) * D_w_inv;

  % PRE-STORED from Phase 2 (no recomputation)
  int_semi = int_semi_store{ts};
  int_semi_q = int_semi_q_store{ts};
  int_semi_w = D_w * int_semi * D_w_inv;
  int_semi_w_q = D_w * int_semi_q * D_w_inv;

  psi_vec_ts = psi_vec_store{ts};
  c_out_per_mode_ts = c_out_store{ts};

  % Coupling vectors (OK: uses q-versions)
  c_in_w = D_w * (abs(int_semi_q) * (3 * psi_vec_ts));
  VinfqJ_ts = VinfqJ_vec(ts);
  VinfJ_ts = VinfJ_vec(ts);
  c_out_w = (D_w_inv * (VinfqJ_ts * 3 * c_out_per_mode_ts))';

  g22_ts = g22_store(ts);

  % Schur complement
  schur = intval(1) - g22_ts - c_out_w * c_in_w;
  if ~(inf(schur) > 0)
    schur_ok = false; break;
  end
  s_inv = 1 / schur;

  % B_dir = (I - G_dir)^{-1}
  B_dir = intval(zeros(m1m2+1));
  B_dir(1:m1m2, 1:m1m2) = intval(eye(m1m2)) + c_in_w * (s_inv * c_out_w);
  B_dir(1:m1m2, m1m2+1) = c_in_w * s_inv;
  B_dir(m1m2+1, 1:m1m2) = s_inv * c_out_w;
  B_dir(m1m2+1, m1m2+1) = s_inv;

  % Base matrices
  base_sup = intval(zeros(m1m2+1));
  base_sup(1:m1m2, 1:m1m2) = Phi_sup_w;
  base_sup(m1m2+1, m1m2+1) = intval(1);

  base_vloc = intval(zeros(m1m2+1));
  base_vloc(1:m1m2, 1:m1m2) = int_semi_w;
  base_vloc(m1m2+1, m1m2+1) = VinfJ_ts;

  base_vloc_q = intval(zeros(m1m2+1));
  base_vloc_q(1:m1m2, 1:m1m2) = int_semi_w_q;
  base_vloc_q(m1m2+1, m1m2+1) = VinfqJ_ts;

  sup_res{ii} = B_dir * base_sup;
  vloc_res{ii} = B_dir * base_vloc;
  vloc_res_q{ii} = B_dir * base_vloc_q;

  Phi_endpt_w_store{ii} = Phi_endpt_w;
  sigma_store(ii) = sigma_step;
  % The endpoint remainder U_ts-E_ts is bounded from the WITHIN-step
  % observation, not from a relative correction multiplied by E_ts.
  G_dir=intval(zeros(m1m2+1));G_dir(1:m1m2,end)=c_in_w;
  G_dir(end,1:m1m2)=c_out_w;G_dir(end,end)=g22_ts;
  remainder_blocks(:,:,ii)=getting_block_norms(G_dir*sup_res{ii},m1m2);
end

if ~schur_ok
  ok = false; seg_r = 0; seg_Z = Inf; Y_max = Inf; U11_K1 = Inf; r_end = 0;
  return;
end

%% Step 2: Signed bare products + additive two-block corrections
[correction,~,~]=getting_multistep_correction(Phi_endpt_w_store,sigma_store,remainder_blocks);
U_11=intval(zeros(K));U_12=U_11;U_21=U_11;U_22=U_11;
V_11=U_11;V_12=U_11;V_21=U_11;V_22=U_11;
Vq_11=U_11;Vq_12=U_11;Vq_21=U_11;Vq_22=U_11;
for ii=1:K
  if strcmp(getenv('GE_DIAGNOSTICS'),'1') && mod(ii,20)==0,fprintf('phase4 composition %d/%d\n',ii,K);end
  S2=getting_block_norms(sup_res{ii},m1m2);
  v2=getting_block_norms(vloc_res{ii},m1m2);
  q2=getting_block_norms(vloc_res_q{ii},m1m2);
  U_11(ii,ii)=S2(1,1);U_12(ii,ii)=S2(1,2);U_21(ii,ii)=S2(2,1);U_22(ii,ii)=S2(2,2);
  V_11(ii,ii)=v2(1,1);V_12(ii,ii)=v2(1,2);V_21(ii,ii)=v2(2,1);V_22(ii,ii)=v2(2,2);
  Vq_11(ii,ii)=q2(1,1);Vq_12(ii,ii)=q2(1,2);Vq_21(ii,ii)=q2(2,1);Vq_22(ii,ii)=q2(2,2);
  Phi_multi=intval(eye(m1m2));sigma_multi=intval(1);
  for j=ii-1:-1:1
    R_bare=intval(zeros(m1m2+1));R_bare(1:m1m2,1:m1m2)=intval(mag(Phi_multi));
    R_bare(end,end)=abs(sigma_multi);
    if j==ii-1,Delta=intval(zeros(2));
    else,Delta=correction{ii-1}(:,2*j+1:2*j+2);end
    v2=getting_block_norms(sup_res{ii}*R_bare*vloc_res{j},m1m2) ...
      +S2*Delta*getting_block_norms(vloc_res{j},m1m2);
    q2=getting_block_norms(sup_res{ii}*R_bare*vloc_res_q{j},m1m2) ...
      +S2*Delta*getting_block_norms(vloc_res_q{j},m1m2);
    V_11(ii,j)=v2(1,1);V_12(ii,j)=v2(1,2);V_21(ii,j)=v2(2,1);V_22(ii,j)=v2(2,2);
    Vq_11(ii,j)=q2(1,1);Vq_12(ii,j)=q2(1,2);Vq_21(ii,j)=q2(2,1);Vq_22(ii,j)=q2(2,2);
    Phi_multi=Phi_multi*Phi_endpt_w_store{j};sigma_multi=sigma_multi*sigma_store(j);
    R_bare(1:m1m2,1:m1m2)=intval(mag(Phi_multi));R_bare(end,end)=abs(sigma_multi);
    u2=getting_block_norms(sup_res{ii}*R_bare,m1m2) ...
      +S2*correction{ii-1}(:,2*j-1:2*j);
    U_11(ii,j)=u2(1,1);U_12(ii,j)=u2(1,2);U_21(ii,j)=u2(2,1);U_22(ii,j)=u2(2,2);
  end
end
U11_K1=U_11(K,1);
verification.composition_verified=true; % Conditional on the supplied one-step enclosures.

%% Step 3: Y and W_Yq (intval throughout)
Y_sup = intval(zeros(K,1));
WYq_sup = intval(zeros(K));

for ii = 1:K
  ts = K_start + ii - 1;
  Yts = intval(0);
  for j = 1:ii
    ts_j = K_start + j - 1;
    F2f = F2ba_fin_vec(ts_j); F2i = F2ba_inf_vec(ts_j);
    F1f = F1ba_fin_vec(ts_j); F1i = F1ba_inf_vec(ts_j);

    Yts = Yts ...
      + U_11(ii,j)*F2f + U_12(ii,j)*F2i ...
      + U_21(ii,j)*F2f + U_22(ii,j)*F2i ...
      + V_11(ii,j)*F1f + V_12(ii,j)*F1i ...
      + V_21(ii,j)*F1f + V_22(ii,j)*F1i;

    WYq_sup(ii,j) = max(Vq_11(ii,j)+Vq_21(ii,j), Vq_12(ii,j)+Vq_22(ii,j));
  end
  if r_prev > 0
    Yts = Yts + (U_11(ii,1) + U_12(ii,1) + U_21(ii,1) + U_22(ii,1)) * r_prev;
  end
  Y_sup(ii) = Yts;
end
Y_max = max(sup(Y_sup));

%% Step 4: Triangular cubic candidates, then the full product-space tests
ban_sup_iv=intval(banorm_vec(K_start:K_end));r_iv=intval(zeros(K,1));
ok=false;seg_r=0;seg_Z=Inf;r_end=0;
for ii=1:K
  past=intval(0);
  if ii>1
    past=WYq_sup(ii,1:ii-1)*((3*ban_sup_iv(1:ii-1)+r_iv(1:ii-1)).*r_iv(1:ii-1).^2);
  end
  local=test_cubic_radius(Y_sup(ii)+past,WYq_sup(ii,ii),ban_sup_iv(ii));
  if ~local.passed,return;end
  r_iv(ii)=intval(local.radius);
end
Pr=Y_sup+WYq_sup*((3*ban_sup_iv+r_iv).*r_iv.^2)-r_iv;
Z_fn=WYq_sup*(3*(2*ban_sup_iv+r_iv).*r_iv);
if all(sup(Pr)<0) && all(sup(Z_fn)<1)
  ok=true;seg_r=max(sup(r_iv));seg_Z=max(sup(Z_fn));r_end=sup(r_iv(end));
  verification.radii_verified=true;
end
end
