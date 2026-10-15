function P4 = pbvp_phase4(bUn_vec, bUn_endpt_vec, r1_vec, ...
    int_semi_store, int_semi_q_store, psi_vec_store, c_out_store, g22_store, ...
    VinfJ_vec, VinfqJ_vec, F1ba_fin_vec, F1ba_inf_vec, F2ba_fin_vec, F2ba_inf_vec, ...
    banorm_vec, tau_vec, mu_m, m, nu)
% Phase 4 for the projected BVP: the same composition as 2D-OK/try_phase4.m
% (signed bare products + additive endpoint remainders U_k-E_k bounded from
% the within-step observation; no relative correction (U-E)E^{-1}).
% Unlike try_phase4, all block tables are returned, since the BVP needs
% U and V for the linearized flows and every row of W_Yq for the X-norm.

m1 = m(1); m2 = m(2); m1m2 = m1*m2;
K = numel(bUn_vec);

[kx_w, ky_w] = ndgrid(0:m1-1, 0:m2-1);
alp_2d = 4*ones(m1, m2); alp_2d(1,:) = 2; alp_2d(:,1) = 2; alp_2d(1,1) = 1;
w_fin = alp_2d(:) .* intval(nu).^kx_w(:) .* intval(nu).^ky_w(:);
D_w = diag(intval(w_fin));
D_w_inv = diag(1 ./ intval(w_fin));

%% Step 1: local directional resolvents and endpoint remainders
sup_res = cell(K,1); vloc_res = cell(K,1); vloc_res_q = cell(K,1);
Phi_endpt_w_store = cell(K,1);
sigma_store = intval(zeros(K,1));
remainder_blocks = intval(zeros(2,2,K));
schur_store = intval(zeros(K,1));
for ts = 1:K
  sigma_step = exp(intval(mu_m) * intval(tau_vec(ts)));
  bUn = bUn_vec{ts};
  ps = intval(zeros(m1m2, m1m2));
  for jj = 1:m1m2
    ps(:,jj) = chebmag_oneside(bUn(:,:,jj));
  end
  Phi_sup_w = D_w * (ps + r1_vec(ts,:)) * D_w_inv;
  Phi_endpt_w = D_w * intval(bUn_endpt_vec{ts}) * D_w_inv;

  int_semi_w = D_w * int_semi_store{ts} * D_w_inv;
  int_semi_w_q = D_w * int_semi_q_store{ts} * D_w_inv;
  c_in_w = D_w * (abs(int_semi_q_store{ts}) * (3 * psi_vec_store{ts}));
  VinfqJ_ts = VinfqJ_vec(ts); VinfJ_ts = VinfJ_vec(ts);
  c_out_w = (D_w_inv * (VinfqJ_ts * 3 * c_out_store{ts}))';
  g22_ts = g22_store(ts);

  schur = intval(1) - g22_ts - c_out_w * c_in_w;
  schur_store(ts) = schur;
  if ~(inf(schur) > 0)
    error('pbvp_phase4:Schur','Schur complement is not strictly positive at step %d.', ts);
  end
  s_inv = 1 / schur;
  B_dir = intval(zeros(m1m2+1));
  B_dir(1:m1m2, 1:m1m2) = intval(eye(m1m2)) + c_in_w * (s_inv * c_out_w);
  B_dir(1:m1m2, m1m2+1) = c_in_w * s_inv;
  B_dir(m1m2+1, 1:m1m2) = s_inv * c_out_w;
  B_dir(m1m2+1, m1m2+1) = s_inv;

  base_sup = intval(zeros(m1m2+1));
  base_sup(1:m1m2, 1:m1m2) = Phi_sup_w; base_sup(end,end) = intval(1);
  base_vloc = intval(zeros(m1m2+1));
  base_vloc(1:m1m2, 1:m1m2) = int_semi_w; base_vloc(end,end) = VinfJ_ts;
  base_vloc_q = intval(zeros(m1m2+1));
  base_vloc_q(1:m1m2, 1:m1m2) = int_semi_w_q; base_vloc_q(end,end) = VinfqJ_ts;

  sup_res{ts} = B_dir * base_sup;
  vloc_res{ts} = B_dir * base_vloc;
  vloc_res_q{ts} = B_dir * base_vloc_q;
  Phi_endpt_w_store{ts} = Phi_endpt_w;
  sigma_store(ts) = sigma_step;
  % The endpoint remainder U_ts-E_ts is bounded from the WITHIN-step
  % observation, not from a relative correction multiplied by E_ts.
  G_dir = intval(zeros(m1m2+1)); G_dir(1:m1m2,end) = c_in_w;
  G_dir(end,1:m1m2) = c_out_w; G_dir(end,end) = g22_ts;
  remainder_blocks(:,:,ts) = getting_block_norms(G_dir*sup_res{ts}, m1m2);
end

%% Step 2: signed bare products + additive two-block corrections
[correction,~,~] = getting_multistep_correction(Phi_endpt_w_store, sigma_store, remainder_blocks);
U_11 = intval(zeros(K)); U_12 = U_11; U_21 = U_11; U_22 = U_11;
V_11 = U_11; V_12 = U_11; V_21 = U_11; V_22 = U_11;
Vq_11 = U_11; Vq_12 = U_11; Vq_21 = U_11; Vq_22 = U_11;
for ii = 1:K
  S2 = getting_block_norms(sup_res{ii}, m1m2);
  v2 = getting_block_norms(vloc_res{ii}, m1m2);
  q2 = getting_block_norms(vloc_res_q{ii}, m1m2);
  U_11(ii,ii)=S2(1,1);U_12(ii,ii)=S2(1,2);U_21(ii,ii)=S2(2,1);U_22(ii,ii)=S2(2,2);
  V_11(ii,ii)=v2(1,1);V_12(ii,ii)=v2(1,2);V_21(ii,ii)=v2(2,1);V_22(ii,ii)=v2(2,2);
  Vq_11(ii,ii)=q2(1,1);Vq_12(ii,ii)=q2(1,2);Vq_21(ii,ii)=q2(2,1);Vq_22(ii,ii)=q2(2,2);
  Phi_multi = intval(eye(m1m2)); sigma_multi = intval(1);
  for j = ii-1:-1:1
    R_bare = intval(zeros(m1m2+1)); R_bare(1:m1m2,1:m1m2) = intval(mag(Phi_multi));
    R_bare(end,end) = abs(sigma_multi);
    if j == ii-1, Delta = intval(zeros(2));
    else, Delta = correction{ii-1}(:,2*j+1:2*j+2); end
    v2 = getting_block_norms(sup_res{ii}*R_bare*vloc_res{j}, m1m2) ...
      + S2*Delta*getting_block_norms(vloc_res{j}, m1m2);
    q2 = getting_block_norms(sup_res{ii}*R_bare*vloc_res_q{j}, m1m2) ...
      + S2*Delta*getting_block_norms(vloc_res_q{j}, m1m2);
    V_11(ii,j)=v2(1,1);V_12(ii,j)=v2(1,2);V_21(ii,j)=v2(2,1);V_22(ii,j)=v2(2,2);
    Vq_11(ii,j)=q2(1,1);Vq_12(ii,j)=q2(1,2);Vq_21(ii,j)=q2(2,1);Vq_22(ii,j)=q2(2,2);
    Phi_multi = Phi_multi*Phi_endpt_w_store{j}; sigma_multi = sigma_multi*sigma_store(j);
    R_bare(1:m1m2,1:m1m2) = intval(mag(Phi_multi)); R_bare(end,end) = abs(sigma_multi);
    u2 = getting_block_norms(sup_res{ii}*R_bare, m1m2) ...
      + S2*correction{ii-1}(:,2*j-1:2*j);
    U_11(ii,j)=u2(1,1);U_12(ii,j)=u2(1,2);U_21(ii,j)=u2(2,1);U_22(ii,j)=u2(2,2);
  end
end

%% Step 3: Y_m and the lumped tables
Y_sup = intval(zeros(K,1));
for ii = 1:K
  Yts = intval(0);
  for j = 1:ii
    Yts = Yts ...
      + U_11(ii,j)*F2ba_fin_vec(j) + U_12(ii,j)*F2ba_inf_vec(j) ...
      + U_21(ii,j)*F2ba_fin_vec(j) + U_22(ii,j)*F2ba_inf_vec(j) ...
      + V_11(ii,j)*F1ba_fin_vec(j) + V_12(ii,j)*F1ba_inf_vec(j) ...
      + V_21(ii,j)*F1ba_fin_vec(j) + V_22(ii,j)*F1ba_inf_vec(j);
  end
  Y_sup(ii) = Yts;
end
lower = tril(true(K));
W_ell = max(U_11 + U_21, U_12 + U_22); W_ell(~lower) = 0;
W_Y = max(V_11 + V_21, V_12 + V_22);   W_Y(~lower) = 0;
W_Yq = max(Vq_11 + Vq_21, Vq_12 + Vq_22); W_Yq(~lower) = 0;

%% Step 4: IVP radii (triangular candidates, then the full product test)
ban = intval(banorm_vec(:)); r_iv = intval(zeros(K,1));
ivp = struct('passed',false,'r',NaN(K,1),'max_r',NaN,'max_Z',Inf);
for ii = 1:K
  past = intval(0);
  if ii > 1
    past = W_Yq(ii,1:ii-1)*((3*ban(1:ii-1)+r_iv(1:ii-1)).*r_iv(1:ii-1).^2);
  end
  local = test_cubic_radius(Y_sup(ii)+past, W_Yq(ii,ii), ban(ii));
  if ~local.passed, break; end
  r_iv(ii) = intval(local.radius);
end
if all(sup(r_iv) > 0)
  Pr = Y_sup + W_Yq*((3*ban+r_iv).*r_iv.^2) - r_iv;
  Z_fn = W_Yq*(3*(2*ban+r_iv).*r_iv);
  ivp.passed = all(sup(Pr) < 0) && all(sup(Z_fn) < 1);
  ivp.r = sup(r_iv); ivp.max_r = max(sup(r_iv)); ivp.max_Z = max(sup(Z_fn));
  ivp.Pr_sup = sup(Pr);
end

P4 = struct('U_11',U_11,'U_12',U_12,'U_21',U_21,'U_22',U_22, ...
  'V_11',V_11,'V_12',V_12,'V_21',V_21,'V_22',V_22, ...
  'Vq_11',Vq_11,'Vq_12',Vq_12,'Vq_21',Vq_21,'Vq_22',Vq_22, ...
  'W_ell',W_ell,'W_Y',W_Y,'W_Yq',W_Yq,'Y_sup',Y_sup, ...
  'schur',schur_store,'remainder_blocks',remainder_blocks,'ivp',ivp);
end
