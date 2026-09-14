# R7P-046--056 — certified boundary-Ising closure

Date: 2026-09-14.

## Scope

This result concerns only the exact four-state even-parity boundary-Ising model inherited from R7P-041--045. It does **not** prove the off-face four-amplitude ceiling and does not transfer to the full seven-coordinate Hessian.

## R7P-046: exact compactification and boundary strata

With `X,Y,Z>=1`, put `r=1/X`, `s=1/Y`, `t=1/Z`. After division by the first positive monomial, the four unnormalised weights are exactly

`(1, 2 s t^3, r t^4, 2 r s t)`

on the closed cube `[0,1]^3`. Hence the cube is the exact physical closure, not the weaker polynomial relaxation. The only zero-probability supports are `{1,2}`, `{1,3}`, and `{1}`. Their covariance ranks are at most one, so `lambda2=0<sigma_*`.

## R7P-047: domain negative controls

For each of the three physical inequalities, a rational witness was found that violates only that inequality while retaining the other two. The R7P-044 shifted-characteristic criterion is evaluated with accepted strict spectral intervals and gives `P(sigma_*)>0` and `P'(sigma_*)<0` in all three cases. Therefore each relaxed model really contains `lambda2>sigma_*`; these witnesses are controls, not physical counterexamples.

## R7P-048--051: frozen polynomial criterion and pilot

Let `d=1+2 s t^3+r t^4+2 r s t`. Positive clearing factors define polynomial numerators

`A=(24 lambda3)^3 d^4 P(sigma_*)`,

`B=(24 lambda3)^2 d^3 P'(sigma_*)`.

By R7P-044, the target is the disjunction `A<=0 OR B>=0`. Exact rational interval Bernstein coefficients are generated coefficient-by-coefficient from strict spectral intervals. The first bounded pilot produced only four unresolved boxes, all in the declared double-root neighborhood; no remote unresolved component survived.

## R7P-052: local equality-neighborhood theorem

At the exact double root `r_*=(1-tau)/(1+tau), s=t=1` with

`tau^2=((2 lambda3-lambda4)(2 lambda3-lambda5))/(4 lambda3^2)`,

exact symbolic reduction proves `A=Ar=As=At=0`. On the dyadic box

`r in [51/128,103/256], s,t in [255/256,1]`,

the interval Hessian satisfies `Arr<0` and the three Schur numerators `Nuu,Nuv,Nvv>0` in tangent coordinates `x=r-r_*`, `u=1-s>=0`, `v=1-t>=0`. Consequently the Hessian quadratic form is strictly negative on every nonzero physical cone direction. Integral Taylor along the convex segment from the double root gives `A<0` away from the root and `A=0` at it. No ordered eigenvalue is differentiated through the degeneracy.

## R7P-053--054: complete cover and independent replay

The refined cover has 436 terminal leaves: 285 strict `A<0`, 150 strict `B>0`, and one leaf discharged by R7P-052. There are zero unresolved leaves. An independent checker, which does not import `boundary_cover.py`, rebuilds the frozen polynomials, interval Bernstein coefficients, the complete split trie, every leaf bound, and the local Hessian certificate. It returns `CHECK_PASS`. Mutation tests reject an altered spectral interval, a corrupted split, a removed leaf, and an invalid leaf-reason substitution.

## R7P-055 theorem

**Certified theorem.** For every probability distribution in the exact four-state boundary-Ising exponential-family closure,

`lambda2(Cov) <= sigma_*`.

Equality is unique and occurs at the exact R7P-045 double-root probability. In compact coordinates it has `s=t=1` and `r=(1-tau)/(1+tau)`. The probabilities are strictly positive inside the four-state boundary model. However, that four-state model itself is the even-parity saturation of the original four-amplitude family, so the equality corresponds to `J6 -> +infinity`, not finite `J6`.

Uniqueness follows because all exterior `A` leaves have strict `A<0`, all `B` leaves have strict `B>0` (which is incompatible with `sigma_*` being the middle eigenvalue), and the local certificate has strict `A<0` except at the exact root.

## R7P-056 dependency impact

The G-lane boundary premise is now paid. H may proceed with R7P-057--064. The off-face I lane is **not** thereby solved: R7P-065 still waits for the intraparity result R7P-064, and the global four-amplitude/full-7D conclusions remain open.
