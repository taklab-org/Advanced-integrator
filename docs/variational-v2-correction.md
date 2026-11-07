# Variational enclosure correction, 2026-11-06

The original v2 finite-section residual omitted the infinite inverse's feedback
from Chebyshev tail rows into finite rows. The corrected Y bound applies the
inverse to the full supported differential residual before projecting. Z0 uses
P K P on a larger window with a certified finite-section remainder. The tail
operator estimates bound output-side rows over all input columns; the Neumann
remainder is included explicitly. See compute_bounds_cheb.m and compute_Y_cheb.m.

PBVP now selects CodeVariational_v2, using the same corrected verifier as GE.
The eigenvalues used in the interval proof are evaluated with INTLAB, including
all arithmetic in epsilon, sigma and L. Newton's floating-point operators are
candidate generators only. Approximate midpoint-based deduplication was removed.
Only the shifted Ohta-Kawasaki cubic at Chebyshev weight one is accepted.

The PDE parameters retain the original exact binary64 input values 0.39 and 1.1;
this revision does not change the IVPs to exact decimal parameters. The PBVP
phase radius remains the explicitly enclosed exact decimal 0.01.

The regression test tests nonlinear coefficient indexing and compares low-order
radii and operator bounds against independently assembled larger-window double
calculations. The latter are diagnostics, not standalone theorem certificates.

Historical outputs and reported radii do not certify the corrected revision.
New PBVP, num5 and num7 computations and the separate fresh-process inspection
have passed; see results.md and validation.md for their scope and results.
All downstream bounds must be recalculated; old Phase 3/4 or intermediate results
are not carried over. The new campaign starts from the public initial data and
regenerates Phase 2 as well, avoiding ambiguity about dataset provenance.
