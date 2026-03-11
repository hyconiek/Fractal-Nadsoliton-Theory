# T152 Current Strict Nad12-Sigma Strict-Source Shannon-Weighted Theta-Export Refinement Candidate Spec

Status: `T152_CURRENT_STRICT_NAD12_SIGMA_STRICT_SOURCE_SHANNON_WEIGHTED_THETA_EXPORT_REFINEMENT_CANDIDATE_SPEC_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

`T81/N346` already show that the repo exports a Shannon-weighted theta-export
refinement candidate only with a canonical-ontology-supported `4 ln 2`
normalizer (and a candidate sigma-int input).

After `T144/N420` and `T149/N418`, the repo now exports strict-side source
upgrades for both needed scalars, and after `T151/N430` it exports the strict
source-upgraded Shannon-weighted feeder-law refinement candidate.

`T152` asks the narrow next question:

```text
can the repo export one strict-source Shannon-weighted theta-export refinement
candidate, i.e. the same refinement intent as T81 but with explicit strict-side
provenance for the Shannon weight and sigma-int input,
without claiming any actual theta export, pair population, feeder support,
residual bridge/export-map object support, or loop break?
```

## Scope

`T152` is scoped only to packaging a refinement candidate above:

1. `N333`
   - actual packaged theta-export candidate,
2. `N430`
   - actual packaged strict-source Shannon-weighted nonequality feeder-law
     refinement candidate,
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

## Candidate theta-export refinement form (strict-source version)

The strongest admissible export form at this stage, if any, is only:

```math
\boldsymbol{\theta}^{cand,sh,strict}_{pair}
:=
\begin{pmatrix}
\theta^{cand,sh,strict}_1 \\
\theta^{cand,sh,strict}_2
\end{pmatrix}
=
\mathcal{M}^{cand,sh,strict}_{nad12,\sigma,res}
\!\left(
\omega,\phi,\sigma^{input}_{int};\,\alpha^{strict}_{geo}
\right)
```

with the exact reading:

1. `\mathcal{M}^{cand,sh,strict}_{nad12,\sigma,res}` is still candidate-only,
2. `\alpha^{strict}_{geo}` is strict-derived:
   `alpha_geo_strict_derived_v1 := 4 ln 2` (`F309/N420`),
3. `\sigma^{input}_{int}` admits the strict-source instantiation:
   `sigma_int_input := sigma_int_strict_derived_v1` (`F307/N418`),
4. `\theta^{cand,sh,strict}_1`, `\theta^{cand,sh,strict}_2` are still
   candidate-only downstream slot values,
5. the pair-indexed target role is already typed by `R1`,
6. the basis-pair population role is already packet-ready only through
   `C48/C49`,
7. no actual `theta_1`, `theta_2` are exported.

## Target

If the answer is positive only at refinement-candidate level, export one strict
source-upgraded packaged refinement candidate:

```text
ThetaPair_strict_nad12_sigma_residual_nonequality_shannon_weighted_export_refinement_candidate_v1
```

with the intended meaning:

```text
actual packaged strict-source Shannon-weighted theta-export refinement candidate
for a future pair-indexed theta supply
on the nad12-sigma residual route

above theta-export-candidate-only language (N333)
above strict-source Shannon feeder-law refinement candidate language (N430)
below actual theta export
below actual pair population
below actual feeder support
below actual residual bridge/export-map object support
below actual loop break
```

## Hard limits

`T152` must not claim:

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

