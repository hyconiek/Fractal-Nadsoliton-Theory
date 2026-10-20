# FIN after KL action — results of one new bounded research atom

**19 October 2026.** Target: a typed, falsifiable minimum history object capable of reconstructing dynamics in one already admitted strict Z12 operator class, without silently introducing an intrinsic source law. Read `PAPER.md` for the full proof.

## Results

- **H1 — CONDITIONAL_THEOREM.** Six labelled normalized stationary Fourier-mode two-time autocorrelations at a common supplied lag determine all six radial jump rates of the strictly positive, reversible D12-symmetric one-label CTMC. Exact 6x6 determinant `-3456 sqrt(3)` and inverse verified symbolically.
- **H2 — EXACT_CLASS_COUNTERFAMILY.** Any omitted mode leaves a positive local alternative. Specifically `w->w+eta*(-1,1,-1,1,-1,1)` preserves complete autocorrelation histories of modes 1–5 at all lags but moves mode 6 by `-12 eta`. At eta=.004 the mode 6 rate changes by -.048.
- **H3 — CONDITIONAL_STABILITY_BOUND.** Max rate error <= `(11/12)*eta/[tau*(Cmin-eta)]` for pointwise correlation noise eta<Cmin. Numerically at tau=.6 and max noise1e-5, actual max rate error 3.1606e-5 against proven bound 6.2287e-5. This is deterministic sensitivity, not sample complexity.
- **H4 — EXACT_FALSIFICATION_CONE.** An alleged six-positive modal correlation vector `C` at one fixed lag belongs to the strictly positive radial Markov class if and only if `B^-1(-log C)>0` componentwise. The synthetic vector `C=(e^-3,e^-1,e^-1,e^-1,e^-1,e^-1)` violates this condition, despite each component being in (0,1).

## Fundamental boundary

The history observables or their model time must be independently sourced. Generating them with a *chosen* CTMC and inferring the same CTMC is an inverse-problem sanity check, **not an origin of dynamics**. Static strict kernel, gradient mobility, GEN-654 KL history energy, A7 count heat bath and the legacy field action remain different typed mathematical objects. No new physical law or source provider is asserted; blocker F01 **OPEN_SOURCE**, F03 **OPEN_SOURCE**, F06 **OPEN_SOURCE**.

## Relation to prior work

GEN-532–554 and 654 proposed least-change/information actions; the immediately preceding work proved their specific Fisher–Rao class obstruction. AGENTS ST8397–8547 already pays several generic graph entropy/refinement and nonuniqueness theorems. Known Markov spectral reconstruction literature overlaps the general inverse method. The local increment consists of the explicit six-mode strict Z12 matrix, invertible/sensitivity certificate, exact positive five-mode counterfamily and log-correlation admissibility inequalities. Do not count these as originating information-time laws.

## Verification

`verify.py` checks all algebra symbolically and reconstructs the correlations with an independent `scipy.linalg.expm` 12x12 matrix exponential over three lags. Numerical values are deterministic synthetic tests from the frozen model, not physical records. Source AGENTS is an uploaded snapshot, not verified to be current live repo. No external file mutations.

**Disposition:** `SIX_MODE_HISTORY_INVERSE_CERTIFICATE_PAID_CONDITIONALLY; FUNDAMENTAL_HISTORY_SOURCE_UNPAID; STOP_FOR_ARCHITECT`.
