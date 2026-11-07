# Corrected results (2026-11-06)

The corrected variational bounds and PBVP Cadiot–Lessard v2 implementation
passed clean recomputation and a separate fresh-process inspection on Windows,
MATLAB R2024b Update 5, INTLAB V13 and Intel MKL (CNR AVX2).
GE used 10 process workers with 2 threads each and a 20-thread client.
All Phase 2/3 data were regenerated from the public initial values.
The decimal values below are rounded reference values; acceptance requires the
strict interval inequalities, rather than bitwise equality between machines.

| Quantity | num5 | num7 |
|---|---:|---:|
| Fourier truncation | [20,20] | [28,28] |
| finite variational block | [7,7] | [7,7] |
| finite endpoint / steps | 13 / 208 | 8 / 128 |
| finite endpoint radius (approx.) | 5.2583e-6 | 2.3792e-6 |
| middle block / maximum step | no middle steps | [12,12] / 0.1 |
| middle steps / terminal entry | 0 / 13 | 69 / 14.9 |
| terminal radius | 0.0018043751029120721 | 0.0060362072956352706 |
| terminal Pr upper | -1.5866289312721146e-7 | -9.6964811499566927e-8 |
| terminal Z upper | 0.12059488597713842 | 0.83938743500267565 |
| complete run time (seconds) | 11004.8205314 | 46187.056853 |

Both GE runs passed all five validation flags. The num7 finite proof ends at
8 directly; it does not reuse or restrict a proof ending at 8.625.

PBVP retains N=[12,12], variational block [3,3], flight time 0.1 and four steps.
Its component radii (theta, phi, u) are
(7.7206292209929278e-12, 3.2453272984168665e-11, 4.8832391016350625e-11),
with maximum Z upper 3.2786313587972981e-9 and positive componentwise margins.
The full-path defect bound is 1.64340258302628e-11; the endpoint bound is
1.5863254085369633e-11. U_DP=4.1963352758018537 uses all step norms.
W_ell(K,1)=1.062299 is a separate table entry and is not identified here with
U_DP or a bound covering the entire trajectory. PBVP run time including
candidate generation was 1113.0605628 seconds.

The independent inspection took 634.8906259 seconds, including 508.5495923
seconds for PBVP revalidation from its saved Newton candidate. It checked GE
source/input hashes, all stored middle attempts and accepted states, histories,
Schur/Neumann/radii conditions, and recomputed the terminal tests. It did not
repeat the expensive GE local-bound generation. PBVP revalidation recomputed
all four steps and confirmed the full-path norm and phase contributions; the
saved radii, Y and Z agreed exactly within this host/run comparison.

Sampled peaks over the GE campaign were 43.65% physical memory load and
18.953 GiB recorded aggregate MATLAB working set (minimum free 72.073 GiB).
These are campaign-wide samples, not per-case memory requirements. Temporary
historical PID reuse was additionally covered by dedicated resource guards.
All computation and inspection processes ended; the three result hashes were
unchanged by inspection.

Result SHA-256 values:

- num5: `a8308ff180b6dfe5a53b96ed43629226d2dc13a3e929539813711db1a72d155d`
- num7: `4397402edf0ec6c8cec79d430d707067d52a85f331ddf1f4d07277a7a9b04ac9`
- PBVP: `6cb8947bec3fd72e26f0d8f13dff3ede76432b45b08168d352fae31f6c7a405e`

Local evidence is retained under outputs/corrected_v2_20261106 (PBVP) and
outputs/corrected_v2_20261106_ge_restart2 (GE and final inspection).
Calculation data and the private INTLAB runtime are excluded from distribution.
The following older results are historical and do not certify this revision.

---

# Recorded results (reference values, not bitwise acceptance targets)

All reference GE results below passed independent final inspection. Reference
interval bounds are summarized here; no binary computation data are included.
Printed decimal values are rounded and are not substitutes for certificates.

| Quantity | num5 | num7, latest recorded tF=8 / middle [12,12] |
|---|---:|---:|
| Fourier truncation | [20,20] | [28,28] |
| finite variational block | [7,7] | [7,7] |
| finite endpoint / used steps | 13 / 208 | 8 / 128 |
| finite endpoint radius | 5.2578e-6 | 2.11658288373e-5 |
| equilibrium norm (approx.) | 0.44104 | 0.73723 |
| full equilibrium defect (approx.) | 4.4325e-8 | 8.2518e-8 |
| middle block / maximum h | no middle steps | [12,12] / 0.1 |
| middle steps / terminal entry time | 0 / 13 | 69 / 14.9 |
| terminal W (approx.) | 1.2017 | 2.1379 |
| terminal V0 / Vq (approx.) | 11.405 / 25.205 | 188.66 / 31.309 |
| terminal radius | 0.00180437436867 | 0.0062164741587029827 |
| terminal Pr upper | -1.58662837435e-7 | -8.42142528481235e-8 |
| terminal Z upper | 0.120594836804 | 0.86456037711573064 |

Reference num5 result SHA-256:
`a7b32762a3abdf8a604883902099d232c34894992bc7395752f4eb83d9285d5a`.
Reference num7 result SHA-256:
`b95022b91aa7dbcacd3b1b46d1adc782a8387c9970a86a2bd2c7011922f38458`.
These identify the original archived results, **not expected hashes for a new run**.

The reference num7 proof first established a full 138-step forward bound on
[0,8.625], then restricted its uniform radius to t=8 and added the polynomial
endpoint discrepancy. The public code instead computes Phase 2/3/4 directly
on [0,8] (128 steps), as requested for this release. All middle bounds are
newly computed from that boundary. The completed corrected direct-eight run
is reported above and is distinct from this historical restriction.

The original num7 rerun took 29019.773 seconds (excluding existing Phase 2/3/4
generation); its independent inspection took 129.466 seconds. Peak physical
memory load was 30% on the 128 GiB host, with MATLAB working set 12.216 GiB.
These are **not total clean reproduction times** or portable memory guarantees.

## PBVP

Recorded serial corrected run: 2026-10-09, Morse index 4 to 3, N=[12,12],
variational block [3,3], time of flight 0.1, phase radius 0.01, four steps.

- IVP radius about 1.643836e-11; IVP Z about 3.442984e-12.
- Component radii (theta, phi, u): approximately
  (8.2301e-12, 3.3385e-11, 5.0973e-11).
- Maximum PBVP Z about 3.4941e-9.
- Strict componentwise margins were positive; run log reported
  (8.2301e-18, 3.3384e-17, 5.0973e-17).
- Recorded runtime: 575.4 seconds, excluding creation of the candidate input.

The PBVP theorem here is a projected BVP result with the specified linear
boundary maps and phase condition. It is not an additional validation of
nonlinear local invariant-manifold parameterizations.

For a new run, regenerated eigenbases and Newton candidates may alter these
numbers considerably. The componentwise inequalities in `run_case` and the
proof script are the acceptance conditions.


## 2026-11-06 variational correction

See [variational-v2-correction.md](variational-v2-correction.md). Earlier numerical
results are historical; the corrected results and inspection are recorded above.
