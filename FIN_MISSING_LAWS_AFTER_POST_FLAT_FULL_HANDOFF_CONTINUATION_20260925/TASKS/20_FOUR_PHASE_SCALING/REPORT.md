# GEOMETRY-FREE-ENERGY-05 — four-phase finite-size scaling

Status: **RECOMPUTED_NUMERIC_SCALING**

## Main result
For the symmetric phase subset `{0,3,6,9}` with n equal replicas of each phase, the hierarchy partition function collapses well in the scaled variable `alpha=beta*n`.  The 95% minimizer concentration threshold scales approximately as `beta95 ~ 4.85/n`.  The heat-capacity-like peak grows rapidly, approximately `n^2` over n=2,4,8, consistent with a first-order-like finite-size signature.

## Core formulas
\[
lpha=eta n,\qquad C_H=eta^2\operatorname{Var}(E).
\]

## Evidence / reproduction
`REPLAYS/replay_four_phase_partition.py` recomputes the minimizer scaling from strict localized-orbit vectors.  Peak values are finite numerical scouts from the same exact recursion.

## Caveats
Only a controlled four-phase subsystem; first-order-like finite-size behavior is not yet a thermodynamic theorem for the full 12-phase alphabet.

## Next question
Extend the same analysis to 8 and 12 phases.
