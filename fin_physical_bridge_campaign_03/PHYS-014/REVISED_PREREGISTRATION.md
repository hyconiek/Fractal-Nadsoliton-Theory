# PHYS-014 — TWO-Q12-HARDWARE-PREREGISTRATION-REVISION

Date: 2026-09-29

Status: **DESIGN_ONLY / VALIDATION UNOPENED**

This supersedes the single-node shortcut in PHYS-010 while preserving the
PHYS-007 ideal 12-bin fingerprint exactly.

## Frozen physical topology

- VRSPAD-144 or an equivalent independently calibrated CTMC platform;
- exactly two contiguous q=12 Potts nodes;
- 24 stochastic channels;
- no other active node;
- g=3 primary condition;
- g=0 negative control;
- observable only: d=(j-i) mod 12.

## Calibration-only channel selection

Before validation, measure per-channel baseline rate and computational
temperature without A7 validation records.

A platform may enter validation only if at least 24 channels satisfy the
predeclared quality box. The current public-data design box is:

    |T - 27.25| <= 1.5
    |baseline residual| <= 1%

The exact selected channel IDs, node split, temperatures, baseline corrections
and all 144 pair weights must be frozen before the first validation record.

## Weight compilation

For each destination channel a:

    E_ab = -(g/2) A_ab
    w_ab = round(T_a E_ab)

Raw stored rows may therefore be asymmetric when T_a differs, while the
intended dimensionless interaction is symmetric.

Do not tune weights against the final 12-bin validation histogram.

## Pre-validation gates

1. `g=0` control including invalid states must pass its frozen distribution
   envelope.
2. Calibration-only CTMC prediction for the exact selected channels and compiled
   weights must be within TV<=0.003 of the ideal FIN fingerprint.
3. Invalid/no-hot/multi-hot events are not discarded; their mass is included
   in the error budget.
4. Drift and serial correlation rules are frozen before validation.
5. No change to A7, g, primary countermodels or the 12-bin observable after
   validation opens.
6. Store raw pair states and timestamps, not only the reduced histogram.

## Validation interpretation

A PASS can establish that a physical stochastic network realizes the
pre-registered N=2 generalized-Curie-Weiss/A7 model to the declared accuracy.

It cannot establish why A7 should occur in nature. That source lane remains
closed unless a new independent source atom is found.

## Public-data warning

PHYS-012 found 3–7% TV discrepancies in the best ordinary published q=4
coupled-Boltzmann output distributions. Therefore no system-level fidelity
claim may be inherited from the publication. The q=12 pair must earn its own
TV<=0.003 calibration gate.
