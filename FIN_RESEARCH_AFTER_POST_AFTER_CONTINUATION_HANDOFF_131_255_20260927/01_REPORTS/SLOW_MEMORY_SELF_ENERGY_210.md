# SLOW-MEMORY-SELF-ENERGY-210
## In the slow Z3 sector, short memory cancels about 80% of the instantaneous projected rate

Date: 2026-09-26

Status:
exact finite-state N=6 symmetry-adapted Mori-Zwanzig calculation.

In the k=4 and k=8 sectors, symmetry makes the problem scalar.

For each slow sector:

    instantaneous projected generator
      A_k
      = -0.113887302674902

    zeroth memory moment
      M0_k
      = +0.090633924632040

    first memory moment
      M1_k
      = 0.028050238847108.

Therefore

    M0_k / |A_k|
      ≈ 0.795821154.

So approximately

    79.582 %

of the naive instantaneous relaxation is canceled by the integrated memory
self-energy.

The first-moment local reduction gives

    lambda_MZ
      =
      (A_k+M0_k)/(1+M1_k)

      = -0.022618912154468.

The exact full microscopic slow eigenvalue from the earlier N=6 calculation is

    lambda_exact
      ≈ -0.022611275187114.

Thus the memory-corrected scalar Z3 mode reproduces the true slow clock to the
previously reported ~0.034%.

## Interpretation

The effective slowness is NOT already present in PLP.

It is generated largely by repeated excursions into unresolved states and
their return.

This is a concrete self-energy mechanism:

    direct projected decay
      +
    hidden recrossing memory
      ->
    much slower effective motion.

Because k=4,8 are pure Z3/base characters under CRT, this renormalization acts
directly on the base coordinate selected by the metastable quotient.
