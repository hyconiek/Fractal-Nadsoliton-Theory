# REF-006 — joint operator/refinement convergence contract

Status: **CONDITIONAL_PASS_WITH_EXACT_NONUNIQUENESS_BOUNDARY**.

## Declared family
On the unit circle use uniform mesh `h=1/M`, measure `h sum_j delta_jh`, and the frozen six-shell operator

`(L_h u)_j = h^-2 sum_{d=1}^6 w_d (2u_j-u_{j+d}-u_{j-d})`

with the strict coefficients from OPERATOR-REFINE-05. No coefficients are refit with M. Let

`kappa = sum_d w_d d^2 = 3.81531411549983`.

For Fourier frequency p the symbol is

`lambda_h(p)=2 h^-2 sum_d w_d(1-cos(p d h))`.

Using `|2(1-cos x)-x^2| <= x^4/12`, for every fixed band `|p|<=P`,

`|lambda_h(p)-kappa p^2| <= h^2 P^4/12 * sum_d w_d d^4`,

where `sum w_d d^4 = 52.41766130594285`. Thus the declared local family has uniform `O(h^2)` symbol convergence on every fixed frequency band. On that band, for `Re z>=rho>0`, the resolvent identity gives a corresponding `O(h^2/rho^2)` response bound. This is enough for a controlled band-limited continuum response and is consistent with the previously observed ~M^-1.98 error slope.

## Boundary / no-go
This does **not** identify the full strict continuum. Two old exact obstructions remain live:

1. visible six-shell data admit both a finite-range continuation with small-k order 2 and a positive `d^-1.8` tail with order 0.8;
2. exact coarse intertwining permits arbitrary nonnegative fiber parameters `mu_n` while leaving the nested coarse sector unchanged.

Therefore metric convergence plus coarse response cannot determine tail order, hidden fiber gaps, temporal memory or the unique continuum law. The fixed-stencil law is a conditional refinement principle, not something selected by FIN.

## Held-out checks
A fine mode with `p h = pi/2` has limiting symbol ratio `0.210703994906` rather than 1, so the fixed-band theorem cannot be promoted to uniform convergence through the Nyquist scale. On mixed nonuniform subdivisions the phrase “offset d=1..6” has no canonical physical-length meaning; a separate irregular refinement law is required. By contrast the nearest-neighbour `c(ell)=kappa0/ell` energy remains subdivision-compatible, showing exactly why the metric backbone transfers more strongly than the full six-shell operator.
