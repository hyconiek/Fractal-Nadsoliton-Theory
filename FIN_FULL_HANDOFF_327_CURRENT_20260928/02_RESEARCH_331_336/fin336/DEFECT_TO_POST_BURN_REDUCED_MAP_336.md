# FIN 336 — DEFECT-TO-POST-BURN REDUCED MAP
## A finite response library replaces repeated full-state preparation propagation

Date: 2026-09-28

Status: **PASS for N=7..10 at the declared burn time.**

## 1. Motivation

Task 331 showed that the exact `D<=6` preparation removes most of the empirical preparation-transfer error, but it still propagated the preparation through the full microscopic burn dynamics.

Task 336 compresses that burn stage itself.

## 2. Construction

For each retained microscopic defect state `x` in the declared J=0 basin with `D<=6`, define the post-burn response vector

`R_x(a) = P_x[ macro outcome at t_burn=24 is a ]`,

where `a` runs over the twelve localized labels plus one explicit `unlocalized` outcome.

Because Markov propagation is linear, for **any** preparation weights `w_x`,

`p_burn(a) = sum_x w_x R_x(a)`.

Thus the dependence on kappa/theta is entirely in the exact defect weights; the expensive microscopic propagation is absorbed once into a small response library.

The responses are computed from the same reversible spectral representations already validated in tasks 311,312,315,317.

## 3. Library sizes

| N | retained J=0 states with D<=6 | full count-space size |
|---:|---:|---:|
| 7 | 2,596 | 31,824 |
| 8 | 5,648 | 75,582 |
| 9 | 9,351 | 167,960 |
| 10 | 11,552 | 352,716 |

Each retained state stores only 13 response probabilities.

The response rows sum to one to better than `1e-12` in all four cases.

## 4. Error structure

For the pinning family, the only preparation truncation error is the omitted tail

`P(D>6 | J=0)`.

Task 328 proves this tail is maximal at kappa=0. Therefore the response library works for every positive pinning strength with the same previously certified tail bound.

At the frozen controller value `theta=2`, equivalently `kappa=2N`, the tails are much smaller:

- N=7: `5.13e-11`;
- N=8: `1.89e-10`;
- N=9: `9.06e-10`;
- N=10: `1.36e-9`.

So for the frozen task333 controller, the finite-defect truncation is practically negligible over N=7..10.

The remaining numerical error is the already-audited spectral burn truncation, not a new fitted preparation residual.

## 5. Scientific meaning

The preparation map now factorizes as

`control parameter -> exact defect weights -> fixed response library -> post-burn prior`.

There is no PCA correction and no separate fit for each kappa.

The map remains N-specific because the response library is generated from the N-specific microscopic generator. A true cross-N analytic response law is still open.

## 6. Verdict

**PASS:** repeated full-state preparation propagation is no longer required for N=7..10 once the response library has been built.

**OPEN:** derive the response vectors analytically/asymptotically from defect geometry so that the library itself need not be recomputed for each N or g.
