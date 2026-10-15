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
newly computed from that boundary. Do not claim its unexecuted full-run output
equals the historical 69-step result.

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
