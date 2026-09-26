# SYMMETRY-REDUCED-LARGEN-CAPACITY-233
## Exact D4 quotient extends deep-core capacity to N=11 and the local exponent moves strongly toward the B4 communication barrier

Date: 2026-09-26

Repository baseline:
main HEAD `ad15a9098ecc5e1282f964ea8b159a8ec608d7c5`.

Status:
exact finite-state potential-theory calculation for the declared leave-one-out
Gibbs chain, using an exact symmetry quotient.

## 1. Symmetry reduction

The capacity problem

    deep Z3 core 0
      vs
    deep Z3 cores 1 union 2

is invariant under the order-8 subgroup

    j -> ±j + 3r,
    r=0,1,2,3.

This subgroup preserves the source sector and swaps the two target sectors.

Therefore the equilibrium potential is constant on its state-space orbits.

The full Dirichlet problem can be solved exactly on the quotient chain.

## 2. State-space reduction

At N=9:

    full count states:
      167,960

    quotient orbits:
      21,100

    reduction:
      7.96x.

At N=10:

    352,716
      ->
    44,352.

At N=11:

    705,432
      ->
    88,410.

The quotient generator remains reversible.

## 3. New exact deep-core capacities

Using the four all-copies-on-one-label deep seeds in each Z3 sector:

    N=9:
      capacity
        = 0.000825425513381

      3*capacity
        = 0.00247627654014

    N=10:
      capacity
        = 0.000450666050755

      3*capacity
        = 0.00135199815227

    N=11:
      capacity
        = 0.000242214388586

      3*capacity
        = 0.000726643165759.

The factor 3 uses the asymptotic/symmetry valley mass 1/3 and is used only for
rate-scale comparison.

## 4. Local exponent

Define

    beta_N
      =
      -log[
        r_(N+1)/r_N
      ].

The latest values are:


    6->7:
      beta_N=0.540297063

    7->8:
      beta_N=0.563392542

    8->9:
      beta_N=0.585833843

    9->10:
      beta_N=0.605172426

    10->11:
      beta_N=0.620903364


The mapped inter-Z3 communication barrier is

    B4=0.662219137127.

Thus the sequence

    0.5634,
    0.5858,
    0.6052,
    0.6209

for N=7->11 moves clearly toward B4 from below.

This is qualitatively different from the earlier N<=8 ambiguity.

## 5. Updated fits

Over N=3,...,11, fixing the exponent to B4 in

    r_N
      =
      A N^alpha exp(-B4 N)

gives

    alpha≈0.628756

with log-RMSE

    0.031326.

If B is fitted freely together with alpha:

    B≈0.628480,
    alpha≈0.420651,

with log-RMSE

    0.027074.

The free exponent has moved substantially upward compared with the earlier
N<=8 fit and is now much closer to B4.

## 6. Boundary

This is strong numerical evidence for an approach toward the B4 communication
scale.

It is NOT yet a proof that

    lim beta_N = B4.

The newest local exponents still show curvature, and simple one-prefactor
extrapolations are not stable.

A proof still requires:
- global communication-path exclusion below B4;
- capacity asymptotics with controlled remainder.

But the finite-N evidence for B4 is now materially stronger.
