# Validation record

The corrected full runs and independent inspection are recorded in the final
section below. The packaging observations immediately below are historical.

## Validation performed while packaging (2026-10-15)

Host: Windows, MATLAB R2024b Update 5, INTLAB V13, Intel MKL with CNR AVX2.
These observations do not claim that another host or BLAS has been tested.

## Completed

- MATLAB Code Analyzer parsing of the selected numerical source trees and
  portable entry/runtime code: passed.
- MATLAB dependency analysis for both GE engines: no dependencies on the
  original project directories. The external numerical dependency is INTLAB;
  MATLAB products detected were MATLAB and Parallel Computing Toolbox.
- Both fixed initial-value constructors reproduce the original `params`,
  `L`, `N` and initial Fourier vector exactly (binary64 input round-trip).
- Fresh GE clients and every member of two-worker process pools passed
  scalar rounding, prepared-cache checks and 512/2352 matrix-enclosure probes.
- Both GE entries started from their code-defined initial values and built
  and Newton-refined the first trajectory candidate at N=20 / N=28.
- **PBVP clean reproduction completed without the original pathdata MAT file**:
  target approximation, stable eigenspace, boundary coordinates, ODE candidate,
  all four rigorous trajectory steps, variational bounds and the final
  componentwise radii inequalities were regenerated.

PBVP reproduction output:

| Quantity | New value |
|---|---:|
| steps | 4 |
| theta radius | 7.5419841017258191e-12 |
| phi radius | 3.1676946751436358e-11 |
| trajectory radius | 4.7613378419379752e-11 |
| maximum Z upper | 3.2026981958393005e-9 |
| minimum component margin | 7.541976487102123e-18 |
| elapsed time, including candidate generation | 789.1768688 seconds |

These differ from the archived reference values and still pass the required
strict inequalities. This is the intended numerical reproducibility policy.

## Not claimed

The GE full-interval calculations were **not rerun to completion during this
packaging task**. Their original archived full proofs and independent final
inspections are the sources for the reference table. The public num7 *direct*
128-step finite proof is a documented orchestration change from the archived
full-138 restriction; no new full-run radii or exact middle count are claimed.

Single-interval GE diagnostic runs were bounded development checks: startup,
rounding, input construction and candidate refinement completed, but the
initial near-origin interval-bound computation was still running when those
checks were ended. They are **not reported as successful Phase2/3 proofs**.
Their private outputs are kept outside the repository and are not distributed.
No mathematical rejection was observed in those bounded checks.

The standard-installation serialized-full-startup worker fallback is provided
for installations lacking the local minimum initializer. This host used its
reviewed minimum initializer; the fallback and other operating systems have
not been separately exercised here. All paths retain the same runtime gates
and must fail if the actual host's scalar or matrix rounding is unsuitable.

After removal of development-only files, the source-only packaging audit,
final MATLAB parser check, source-manifest round-trip and Git whitespace
check all passed. No retained GE dependency was removed by that cleanup.
The source inventory is in `release-files.json`. No full original datasets,
runtime cache, archived process logs, or proprietary toolbox binaries are
part of that inventory.

## MATLAB launcher update (2026-10-16)

The MATLAB-only launcher was exercised on Windows with MATLAB R2024b:
three-example dry-run planning; early rejection of occupied output paths;
paths containing spaces; actual fresh MATLAB child startup; propagation of
a deliberately invalid INTLAB installation as a nonzero child failure; and
unchanged environment/path/working directory in the parent session.

The successful three-child sequence was tested using a private numerical
entry-point fixture. Each child returned a distinct process ID and a fixture
summary in num5/num7/PBVP order. These fixture success markers are test outputs,
not mathematical proof certificates. The fixture is created in a temporary
directory and removed after the tests; it is never used by the public runner.
The launcher tests do not rerun the numerical proofs or extend the numerical
validation claims above. The release audit is now `tests/check_package.m`.


## 2026-11-06 variational correction

See [variational-v2-correction.md](variational-v2-correction.md). Earlier numerical
results are historical; complete corrected-revision validation is recorded below.


## Corrected revision: completed full runs and independent inspection

PBVP, num5 and num7 all passed on 2026-11-06; see results.md for numbers,
run scope, timings and hashes. Low-order enlarged-window diagnostics and the
fresh two-worker initialization regression passed before production. The
low-order double diagnostics are not substitutes for the interval proofs.

GE regenerated Phase 2 and Phase 3, passed the source/data audit and all five
final flags, and checked strict terminal Pr<0 and Z<1. A fresh 20-thread client
passed scalar rounding and 512/2352 matrix-enclosure checks, then reread all GE
proof records and recomputed the saved middle-step algebra and terminal tests.
Every accepted state and its history matched; all local Schur and Neumann
conditions passed. This inspection checks stored local evidence and downstream
algebra; it does not independently regenerate every GE local bound.

PBVP was independently revalidated from its saved Newton candidate, including
all four local variational steps. Its radii, Y and Z matched the saved result.
The X-norm checks use every time step, and all three Z components include their
phase terms. The full-path and endpoint defects were checked separately.

All 612 frozen MATLAB sources remained unchanged during the GE campaign and
inspection. The original PBVP source manifest differs only in the three later
startup repairs (ai_worker_init, ge_verified_pool and intlab_worker_start);
its mathematical sources are unchanged. All result SHA-256 values remained
unchanged. Input parameters preserve the original binary64 values.
No computation or inspection process remained active at final readback.
