# T154 Current Strict Nad12-Sigma Strict-Source Shannon-Weighted Theta-Export Support Refinement Candidate Spec

Status: `T154_CURRENT_STRICT_NAD12_SIGMA_STRICT_SOURCE_SHANNON_WEIGHTED_THETA_EXPORT_SUPPORT_REFINEMENT_CANDIDATE_SPEC_NO_FALSE_PASS`  
As of: `2026-03-11`

## Goal

After `N431`, the nad12-sigma residual route already exports:

```text
one actual packaged strict-source Shannon-weighted theta-export refinement candidate
```

After `N432`, the same route also exports:

```text
one actual packaged strict-source Shannon-weighted pair-population refinement candidate
```

The next honest question is narrower:

```text
can the current repo already export one actual packaged
strict-source Shannon-weighted theta-export support refinement candidate
for the same route,
without pretending that any actual theta export,
actual pair population,
actual feeder support,
actual residual bridge/export-map object support,
or actual loop break is already present?
```

`T154` is intentionally weaker than:

```text
actual theta export
actual pair population
actual feeder support
actual residual bridge/export-map object support
actual sandbox loop break
actual E_orient
strict-core selector closure
ToE closure
```

## Scope

`T154` is scoped only to the following already exported lanes:

1. `N431`
   - actual packaged strict-source Shannon-weighted theta-export refinement candidate,
2. `N432`
   - actual packaged strict-source Shannon-weighted pair-population refinement candidate,
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

## Candidate theta-export support-refinement form (strict-source version)

The strongest admissible export form at this stage, if any, is only:

```math
\boldsymbol{\theta}^{cand,sh,sup,strict}_{pair}
:=
\begin{pmatrix}
\theta^{cand,sh,sup,strict}_1 \\
\theta^{cand,sh,sup,strict}_2
\end{pmatrix}
=
\mathcal{M}^{cand,sh,sup,strict}_{nad12,\sigma,res}
\!\left(
\omega,\phi,\sigma^{input}_{int};\,\alpha^{strict}_{geo}
\right)
```

with the exact reading:

1. `\mathcal{M}^{cand,sh,sup,strict}_{nad12,\sigma,res}` is still one
   candidate-only export relation,
2. the theta candidate is now supported by the already exported strict-source
   Shannon-weighted pair-population refinement layer,
3. `\theta^{cand,sh,sup,strict}_1`, `\theta^{cand,sh,sup,strict}_2` are still
   candidate-only downstream slot values,
4. no actual `theta_1`, `theta_2` are exported.

## Exact decision question

Can the current repo honestly export anything stronger than:

```text
actual packaged strict-source Shannon-weighted theta-export refinement candidate
```

and stronger than:

```text
actual packaged strict-source Shannon-weighted pair-population refinement candidate
```

namely:

```text
one actual packaged strict-source Shannon-weighted theta-export support
refinement candidate
for the nad12-sigma residual route,
using the already exported theta-export refinement,
the already exported pair-population refinement,
the already exported object-support witness layers,
and the already exported pair-indexed target-slot and conditional population
schema,
while remaining explicitly below:
- actual theta export,
- actual pair population,
- actual feeder support,
- actual residual bridge/export-map object support,
- actual loop break?
```

## Target

If the answer is positive only at support-refinement-candidate level, export:

```text
ThetaPair_strict_nad12_sigma_residual_nonequality_shannon_weighted_support_refinement_candidate_v1
```

with the intended meaning:

```text
actual packaged candidate-only strict-source Shannon-weighted theta-export support
refinement layer for a future pair-indexed theta supply
on the nad12-sigma residual route

above strict-source Shannon-weighted theta-export-refinement-candidate-only language
above strict-source Shannon-weighted pair-population-refinement-candidate-only language
below actual theta export
below actual pair population
below actual loop break
```

## Hard limits

`T154` must not claim:

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

