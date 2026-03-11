# T153 Current Strict Nad12-Sigma Strict-Source Shannon-Weighted Pair-Population Refinement Candidate Spec

Status: `T153_CURRENT_STRICT_NAD12_SIGMA_STRICT_SOURCE_SHANNON_WEIGHTED_PAIR_POPULATION_REFINEMENT_CANDIDATE_SPEC_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

`T82/N347` already show that the repo exports a Shannon-weighted pair-population
refinement candidate only with a canonical-ontology-supported `4 ln 2`
normalizer.

After `T144/N420` and `T149/N418`, the repo now exports strict-side source
upgrades for both needed scalars, and after `T152/N431` it exports the strict
source-upgraded Shannon-weighted theta-export refinement candidate.

`T153` asks the narrow next question:

```text
can the repo export one strict-source Shannon-weighted pair-population
refinement candidate, i.e. the same refinement intent as T82 but with explicit
strict-side provenance for the Shannon weight (and upstream strict-source theta
refinement), without claiming any actual pair population, theta export, feeder
support, residual bridge/export-map object support, or loop break?
```

## Scope

`T153` is scoped only to packaging a refinement candidate above:

1. `N334`
   - actual packaged pair-population candidate,
2. `N431`
   - actual packaged strict-source Shannon-weighted theta-export refinement
     candidate,
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

## Candidate pair-population refinement form (strict-source version)

The strongest admissible export form at this stage, if any, is only:

```text
populated_instance^{cand,sh,strict}(\theta_1^{cand,sh,strict},\theta_2^{cand,sh,strict}) := {
  theta_1: theta_1^{cand,sh,strict},
  theta_2: theta_2^{cand,sh,strict},
  u_1: cos(theta_1^{cand,sh,strict})c_1 + sin(theta_1^{cand,sh,strict})s_1,
  u_2: cos(theta_2^{cand,sh,strict})c_2 + sin(theta_2^{cand,sh,strict})s_2,
  orientation_slice_candidate: span{u_1,u_2},
  strict_weight: alpha_geo_strict_derived_v1
}
```

with the exact reading:

1. the population is still candidate-only syntax,
2. `alpha_geo_strict_derived_v1 := 4 ln 2` is strict-derived (`F309/N420`),
3. `theta_1^{cand,sh,strict}`, `theta_2^{cand,sh,strict}` are still
   candidate-only downstream slot values (as in `N431`),
4. no actual populated basis-pair instance is exported,
5. no actual `u_1`, `u_2` are exported.

## Target

If the answer is positive only at refinement-candidate level, export one strict
source-upgraded packaged refinement candidate:

```text
BasisPair_strict_nad12_sigma_residual_nonequality_shannon_weighted_population_refinement_candidate_v1
```

with the intended meaning:

```text
actual packaged strict-source Shannon-weighted pair-population refinement candidate
layer for a future pair-indexed populated basis-pair instance
on the nad12-sigma residual route

above pair-population-candidate-only language (N334)
above strict-source Shannon theta-export refinement candidate language (N431)
below actual pair population
below actual theta export
below actual loop break
```

## Hard limits

`T153` must not claim:

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

