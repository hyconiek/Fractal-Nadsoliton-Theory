# P397 Current Actual Strict Nad12-Sigma Strict-Source Shannon-Weighted Pair-Population Support Refinement Candidate Probe

Status: `P397_CURRENT_ACTUAL_STRICT_NAD12_SIGMA_STRICT_SOURCE_SHANNON_WEIGHTED_PAIR_POPULATION_SUPPORT_REFINEMENT_CANDIDATE_PROBE_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

Check whether the current repo honestly exports one actual packaged strict-source
Shannon-weighted pair-population support refinement candidate for the
nad12-sigma residual route.

This probe is the strict-source analogue of `P329` (canonical `4 ln 2`).

## Probe table

| Check | Verdict | Basis |
|---|---|---|
| strict-source Shannon-weighted pair-population refinement already exported? | YES | `N432` |
| strict-source Shannon-weighted theta-export support refinement already exported? | YES | `N433` |
| residual bridge/export-map object-support witness already exported? | YES | `N343` |
| nad12-sigma object-support support witness already exported? | YES | `N344` |
| pair-indexed residual target-slot language already exported? | YES | `R1` |
| conditional populated-instance schema already packet-ready? | YES | `C48/C49` |
| actual residual bridge/export-map object support exported? | NO | `N302` remains in force |
| actual pair population exported? | NO | current route still candidate-only |
| actual loop break exported? | NO | sandbox `N18` remains in force |

## Verdict

The current repo honestly exports:

```text
BasisPair_strict_nad12_sigma_residual_nonequality_shannon_weighted_support_refinement_candidate_v1
```

with the exact meaning:

```text
actual packaged strict-source Shannon-weighted pair-population support refinement candidate

above strict-source Shannon-weighted pair-population refinement candidate
above strict-source Shannon-weighted theta-export support refinement candidate
below actual pair population
below actual theta export
below actual feeder support
below actual residual bridge/export-map object support
below actual loop break
```

## Explicit non-claims

`P397` does **not** claim:

1. actual nad12-sigma feeder support,
2. actual residual bridge/export-map object support,
3. actual sigma asymmetry law,
4. actual `f(sigma_int)` weighting law,
5. actual `lambda_1`, `lambda_2`,
6. actual `u_1`, `u_2`,
7. actual `theta_1`, `theta_2`,
8. actual pair population,
9. actual loop break,
10. actual `E_orient`,
11. admissible `S_sel_int`,
12. strict-core selector closure,
13. ToE closure.

