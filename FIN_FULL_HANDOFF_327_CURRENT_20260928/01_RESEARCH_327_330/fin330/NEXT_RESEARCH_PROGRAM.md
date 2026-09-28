# FIN — NEXT RESEARCH PROGRAM AFTER 330

Date: 2026-09-28

## P0 — 331 DEFECT-BASED-OBSERVED-PROCESS-PREDICTOR

Replace the empirical preparation map used in the 313-style process predictor by the exact few-defect distribution truncated at `D<=6`.

Acceptance tests:

1. preserve the exact basin condition;
2. keep `unlocalized` as an explicit outcome rather than silently postselecting it away;
3. propagate the exact 328 tail bound into the observed-process TV budget;
4. compare preparation and joint-law errors on already opened N=7..12 without refitting per N;
5. if the model is then changed, freeze it before a genuinely new holdout in `g` or preparation protocol, not merely another N in the same calibration lane.

Kill-test: if exact D<=6 preparation does not reduce or explain `eps_prep`, identify which omitted part of the preparation map creates the residual rather than adding another PCA correction.

## P0/P1 — 332 PREPARATION-CONTROLLER-WITHOUT-BASIN-ORACLE

The hard reflecting J=0 wall is operationally strong. Compare:

- unrestricted biased heat-bath + explicit localization outcome;
- soft confining field based on available observables;
- feedback controller with recorded information/control cost.

Measure escape/rejection probability, mixing time, and downstream error.

## P1 — 333 FREEZE A CONTROL PROTOCOL

Before a new experimental-style holdout, explicitly select one control convention:

- fixed kappa;
- fixed theta=kappa/N;
- fixed KL preparation cost.

Current N<=10 evidence favors fixed theta as the most stable few-defect protocol, but this is not yet a physical law. The selection must be justified operationally and frozen before the new test.

## Formal audit lane — 334 FULL MPFR/ARB REPLAY OF 327

Replay the complete six-dimensional and triple-junction covers with a high-precision directed interval backend such as Arb/MPFR, retaining box logs/checkpoints and exact source hashes.

This is the cleanest route to external/repository acceptance of `Gamma=B4` if that proof standard is required.

## Conditional lane — EYRING-KRAMERS PREFACTOR

Do not export an unconditional prefactor until the global barrier certificate is accepted at the desired audit standard. Conditional calculations may proceed, clearly labelled `IF Gamma=B4`.

## P1 — ANALYTIC INITIAL-SLIP COEFFICIENTS

Continue 323/328 by deriving the backward-observable gradients/Hessians needed for the six initial-slip amplitudes. Compare the resulting large-N coefficients against the exact defect representation in their overlap regime.
