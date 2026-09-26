# COMMUTING-TRANSFORMATION-GEOMETRY-245
## Independent commuting information-preserving transformations generate discrete torus geometry and d-dimensional diffusion

Date: 2026-09-26

Status:
exact finite-group theorem plus numerical heat-kernel replay.

Let

    T_1,...,T_d

be commuting bijections of a finite set V.

Assume:
1. T_mu has finite order n_mu;
2. every vertex can be written uniquely as

       T_1^(a_1)...T_d^(a_d) v_0,

   with

       a_mu in Z_(n_mu).

Equivalently, the transformations generate a free transitive action of

    G
      =
    Z_(n_1) × ... × Z_(n_d).

## 1. Geometry

Connect every vertex v to

    T_mu v
    and
    T_mu^(-1) v

for every mu.

Then the resulting graph is exactly

    boxed:
    C_(n_1)
      square
    ...
      square
    C_(n_d),

the Cartesian product of d cycles.

No coordinates need to be inserted after the transformation action is given.

The group exponents ARE the coordinates.

## 2. Laplacian

For edge rate r_mu in direction mu,

    L
      =
      sum_mu
      r_mu(
        2I-T_mu-T_mu^(-1)
      ).

Characters labelled by

    k=(k_1,...,k_d)

have eigenvalues

    boxed:
    lambda(k)
      =
      2 sum_mu
      r_mu[
        1-cos(
          2 pi k_mu/n_mu
        )
      ].

## 3. Fixed local activity budget

If the total per-site event-participation rate is rho and all 2d ports are
equivalent, then

    r_mu=rho/(2d).

Hence

    boxed:
    lambda(k)
      =
      (rho/d)
      sum_mu[
        1-cos(
          2 pi k_mu/n_mu
        )
      ].

For equal large side length L,

    lambda(k)
      approximately
      (2 pi^2 rho/(d L^2))
      |k|^2

for low k.

This is the standard discrete diffusion spectrum in d dimensions.

## 4. Heat-kernel scaling

For large equal tori and times between the local and finite-size scales,

    P_return(t)
      ~
      t^(-d/2).

A numerical replay for L=64 gives median effective spectral dimensions over
the intermediate window t=10...100:


    d=1:
      d_eff≈1.008258

    d=2:
      d_eff≈2.034241

    d=3:
      d_eff≈3.080138


The small excess above integer d is a finite-size/discrete-time-window effect.

## Boundary

The theorem says:

    sourced commuting transformations
      ->
    sourced discrete geometry.

It does NOT tell FIN how many independent transformations exist.

So it converts the spatial-dimension problem into a transformation-rank/source
problem.
