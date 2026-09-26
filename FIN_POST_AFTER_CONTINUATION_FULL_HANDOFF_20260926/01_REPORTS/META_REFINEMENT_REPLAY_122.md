# META-REFINEMENT-REPLAY-122
## Three distinct recursively sourced fiber laws pass two refinement levels exactly

Date: 2026-09-26

A representative two-scale metastable Laplacian was used with normalized

    k3=100,
    k4=1,

so its slow spectral gap is

    delta0=3.

For the recursively sourced rules

    mu=c delta

the following choices were replayed:


### c=0.5

    mu level 0 = 1.5
    gap level 1 = 3
    mu level 1 = 1.5

    level-1 generator intertwining error
      = 0.000e+00

    level-2 generator intertwining error
      = 0.000e+00

    level-1 semigroup error at t=0.37
      = 6.072e-15

    level-2 semigroup error at t=0.37
      = 1.782e-14

### c=1.0

    mu level 0 = 3
    gap level 1 = 3
    mu level 1 = 3

    level-1 generator intertwining error
      = 0.000e+00

    level-2 generator intertwining error
      = 0.000e+00

    level-1 semigroup error at t=0.37
      = 2.260e-14

    level-2 semigroup error at t=0.37
      = 1.568e-14

### c=2.0

    mu level 0 = 6
    gap level 1 = 3
    mu level 1 = 6

    level-1 generator intertwining error
      = 1.137e-13

    level-2 generator intertwining error
      = 1.392e-13

    level-1 semigroup error at t=0.37
      = 5.480e-15

    level-2 semigroup error at t=0.37
      = 1.473e-14

All errors are numerical roundoff.

Thus an unused second refinement level does not discriminate among
c=1/2,1,2.

The test fails for a structural reason, not because of insufficient numerical
precision: every c belongs to an exact analytic family.
