# Packaging changes relative to the recorded source trees

The source hashes in `source-provenance.json` identify the selected original
files. `release-files.json` identifies this release's actual source bytes.
Neither manifest supplies a numerical theorem certificate.
Text line endings and trailing whitespace are normalized for Git distribution.

1. Case-local numerical kernels retain their original formulas and strict
   tests. Unrelated development demos and tests requiring external datasets
   or developer paths are excluded. Original copyright headers are retained.
2. GE's in-function host-specific INTLAB startup is replaced by `run_case`.
   Absolute paths, report directories, process counts and thread counts are
   supplied by the launcher. Phase 2/3 candidates and bounds are generated
   from the explicitly defined initial Fourier coefficients.
3. num5 keeps N=20, final finite time 13 and 208 fixed steps. A failed target
   is not accepted by silently selecting a smaller final finite time.
4. num7 keeps N=28, equilibrium/terminal support 28, finite variational block
   7 and middle block 12. **Its finite target is changed to 8, with 128 steps**.
   Phase 4 tests this newly generated prefix directly. It no longer uses the
   original full-138 restriction or a saved parent's terminal radius.
5. For GE, the successful independent-column `chebmag` implementation is reused in
   process pools. Arithmetic inside each column is unchanged; output columns
   are collected in index order. PBVP retains its own recorded serial
   Chebyshev kernel. Parallel counts are options.
6. PBVP's input MAT-file load is redirected to a fresh output directory. A
   new numerical generator constructs the origin's four mode eigenvectors,
   recomputes the target approximation and stable eigenspace, refines the
   projected boundary coordinates and integrates a fresh candidate ODE.
   No reference orbit, computed eigenspace or saved inverse is distributed.
   The candidate alone does not prove anything; the existing strict PBVP
   integrator and radii inequalities must still pass.
7. A one-interval GE `Smoke` option exercises real Phase 2/3 calculations.
   It is explicitly labelled as a kernel check, never a full proof.
8. The launcher requires an empty output directory and writes separate configuration, log, results and strict-check
   summaries. The post-run GE check recomputes middle/terminal algebraic
   tests from the new local bounds. It does not import the original private
   final-inspection reports or pretend to repeat their entire archival audit.
9. `reproduce.m` launches one or all three examples from MATLAB, each in a
   fresh child process with its own preferences and explicit environment.
   Commands are passed as argument arrays instead of shell-built strings.
   The former Python launcher and Python packaging audit have been replaced
   by MATLAB functions; neither running nor auditing the release needs Python.

The fixed initial values are stored as 17-significant-digit MATLAB constants,
so the original binary64 inputs round-trip. These constants are part of the
problem definition. This precision requirement for the *input* must not be
confused with bitwise matching of computed error enclosures across machines.
