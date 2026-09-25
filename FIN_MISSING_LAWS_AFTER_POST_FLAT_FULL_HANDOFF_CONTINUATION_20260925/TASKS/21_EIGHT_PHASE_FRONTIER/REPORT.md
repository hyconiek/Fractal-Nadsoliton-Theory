# GEOMETRY-TRANSITION-06 — eight-phase extension

Status: **IN_PROGRESS_NUMERICAL_SCOUT**

## Main result
The symmetric 8-phase alphabet complementary to `{0,3,6,9}` has begun to show the same scaled variable `alpha=beta*n`.  Fresh lightweight replay in this build reproduces close n=1/n=2 partition behavior at fixed alpha; an attempted n=4 full recursion exceeded the short execution budget and is therefore not promoted here.  Earlier conversation-level n=2/n=4 collapse notes remain unverified by this package.

## Core formulas
Current retained statement: `alpha=beta*n` remains the leading scaling candidate; no 8-phase transition exponent or thermodynamic conclusion is exported.

## Evidence / reproduction
`REPLAYS/replay_eight_phase_partial.py` computes n=1 and n=2 values at selected alpha.

## Caveats
IN PROGRESS.  Do not cite an n=4 result from this handoff.

## Next question
Implement symmetry-compressed memoization or transfer-matrix reduction for n>=4.
