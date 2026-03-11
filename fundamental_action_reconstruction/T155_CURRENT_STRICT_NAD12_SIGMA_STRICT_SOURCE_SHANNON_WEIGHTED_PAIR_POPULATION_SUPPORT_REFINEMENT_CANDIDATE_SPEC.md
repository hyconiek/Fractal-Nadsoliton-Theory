# T155 Current Strict Nad12-Sigma Strict-Source Shannon-Weighted Pair-Population Support Refinement Candidate Spec

Status: `T155_CURRENT_STRICT_NAD12_SIGMA_STRICT_SOURCE_SHANNON_WEIGHTED_PAIR_POPULATION_SUPPORT_REFINEMENT_CANDIDATE_SPEC_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

After `N432`, the nad12-sigma residual route already exports:

```text
one actual packaged strict-source Shannon-weighted pair-population refinement candidate
```

After `N433`, the same route also exports:

```text
one actual packaged strict-source Shannon-weighted theta-export support refinement candidate
```

The next honest question is narrower:

```text
can the current repo already export one actual packaged
strict-source Shannon-weighted pair-population support refinement candidate
for the same route,
without pretending that any actual pair population,
actual theta export,
actual feeder support,
actual residual bridge/export-map object support,
or actual loop break is already present?
```

`T155` is intentionally weaker than:

```text
actual pair population
actual theta export
actual feeder support
actual residual bridge/export-map object support
actual sandbox loop break
actual E_orient
strict-core selector closure
ToE closure
```

## Scope

`T155` is scoped only to the following already exported lanes:

1. `N432`
   - actual packaged strict-source Shannon-weighted pair-population refinement candidate,
2. `N433`
   - actual packaged strict-source Shannon-weighted theta-export support refinement candidate,
3. `N343`
   - actual residual bridge/export-map object-support witness,
4. `N344`
   - actual nad12-sigma object-support support witness,
5. `R1`
   - pair-indexed residual target-slot language,
6. `C48`
   - packet-ready minimal basis-pair export skeleton,
7. `C49`
   - packet-ready conditional populated-instance schema,
8. `N302`
   - actual residual bridge/export-map object support still frozen,
9. sandbox `N18`
   - theta/population loop still not actually broken.

## Candidate pair-population support-refinement form (strict-source version)

The strongest admissible export form at this stage, if any, is only:

```text
populated_instance^{cand,sh,sup,strict}_{pair} := {
  theta_1: theta_1^{cand,sh,sup,strict},
  theta_2: theta_2^{cand,sh,sup,strict},
  u_1: cos(theta_1^{cand,sh,sup,strict})c_1 + sin(theta_1^{cand,sh,sup,strict})s_1,
  u_2: cos(theta_2^{cand,sh,sup,strict})c_2 + sin(theta_2^{cand,sh,sup,strict})s_2,
  orientation_slice_candidate: span{u_1,u_2},
  strict_weight: alpha_geo_strict_derived_v1
}
```

with the exact reading:

1. the populated instance remains candidate-only,
2. the weight is strict-derived via `alpha_geo_strict_derived_v1 := 4 ln 2` (`N420`),
3. no actual populated basis-pair instance is exported,
4. no actual `u_1`, `u_2` are exported.

## Target

If the answer is positive only at support-refinement-candidate level, export:

```text
BasisPair_strict_nad12_sigma_residual_nonequality_shannon_weighted_support_refinement_candidate_v1
```

with the intended meaning:

```text
actual packaged candidate-only strict-source Shannon-weighted pair-population support
refinement layer for a future pair-indexed populated basis-pair instance
on the nad12-sigma residual route

above strict-source Shannon-weighted pair-population-refinement-candidate-only language
above strict-source Shannon-weighted theta-export-support-refinement-candidate-only language
below actual pair population
below actual theta export
below actual loop break
```

## Hard limits

`T155` must not claim:

1. actual nad12-sigma feeder support,
2. actual residual bridge/export-map object support,
3. actual sigma asymmetry law,
4. actual `f(sigma_int)` weighting law beyond refinement syntax,
5. actual `lambda_1`, `lambda_2`,
6. actual `u_1`, `u_2`,
7. actual `theta_1`, `theta_2`,
8. actual populated basis-pair instance,
9. actual sandbox loop break,
10. actual `E_orient`,
11. admissible `S_sel_int`,
12. strict-core selector closure,
13. `QW-2191` discharge,
14. ToE closure.

