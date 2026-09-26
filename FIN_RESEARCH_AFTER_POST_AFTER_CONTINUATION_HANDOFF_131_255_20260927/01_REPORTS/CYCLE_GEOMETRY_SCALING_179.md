# CYCLE-GEOMETRY-SCALING-179
## The conditional two-port topology has a controlled one-dimensional graph limit

Date: 2026-09-26

Status:
exact graph spectrum/resistance theorem conditional on the C_n incidence of
report 178.

For uniform edge conductance c, the cycle Laplacian has eigenvalues

    boxed:
    lambda_m
      =
      2c[
        1-cos(2 pi m/n)
      ],
    m=0,...,n-1.

The first positive eigenvalue is

    lambda_1
      =
      4c sin^2(pi/n)
      ~
      4 pi^2 c/n^2.

Thus the low graph-frequency modes develop the standard n^-2 cyclic scaling.

## Effective resistance

For two vertices separated by r edges,

    boxed:
    R(r)
      =
      r(n-r)/(n c).

For r << n,

    R(r)
      approximately
      r/c.

So with the earlier static refinement identification

    1/c proportional to edge length,

short-distance resistance is linear in graph path length.

## Spectral-dimension statement

For large n and intermediate diffusion times

    1 << c t << n^2,

the return probability of the cycle heat kernel obeys the one-dimensional
asymptotic

    P_return(t)
      approximately
      1/sqrt(4 pi c t).

Therefore the graph spectral dimension tends to

    d_s=1

in that scaling regime.

This is a genuine graph-theoretic one-dimensional limit.

It is NOT yet physical space:
- the incidence law is conditional;
- the edge length/unit is unsourced;
- the heat generator is not automatically physical propagation.
