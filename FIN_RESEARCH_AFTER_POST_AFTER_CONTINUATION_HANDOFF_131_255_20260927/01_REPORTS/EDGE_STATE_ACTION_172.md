# EDGE-STATE-ACTION-172
## A minimal gauge-equivariant edge action still cannot select bounded-degree sparsity with fixed couplings

Date: 2026-09-26

Status:
mean-field no-go for the smallest joint activation/phase action.

## 1. Minimal node-edge variables

For every candidate pair x<y introduce:

    b_xy in {0,1}
      edge activation,

    a_xy=-a_yx in S1
      independent link phase.

Node phases transform as

    theta_x -> theta_x+alpha_x,

and link phases as

    a_xy
      ->
      a_xy+alpha_y-alpha_x.

Then the edge alignment

    a_xy-(theta_y-theta_x)

and cycle flux

    Phi_xyz
      =
      a_xy+a_yz+a_zx

are gauge invariant.

A minimal permutation/gauge-equivariant energy is schematically

    H
      =
      mu sum_e b_e
      -h sum_(xy)
         b_xy cos[a_xy-(theta_y-theta_x)]
      -J sum_(xyz)
         b_xy b_yz b_zx cos Phi_xyz.

This is the smallest model combining:
- activation;
- relative phase;
- nontrivial holonomy.

## 2. Symmetric mean-field activation equation

Take a homogeneous phase-aligned trial state with

    E[b_xy]=p

for every pair.

There are

    P=n(n-1)/2

candidate edges and

    T=n(n-1)(n-2)/6

candidate triangles.

Absorb the bounded edge-alignment contribution into an O(1) effective
chemical potential mu_eff.

The mean-field free energy is

    F(p)
      =
      P[
        mu_eff p
        +p log p
        +(1-p)log(1-p)
      ]
      -J T p^3.

Stationarity gives

    boxed:
    log[p/(1-p)]
      =
      -mu_eff
      +J(n-2)p^2.

## 3. Test finite mean degree

Bounded mean degree d requires

    p=d/(n-1)
      ~d/n.

Then

    log[p/(1-p)]
      =
      -log n
      +log d
      +o(1),

while for fixed J

    J(n-2)p^2
      =
      O(1/n)
      ->0.

The right-hand side therefore remains O(1) for fixed mu_eff,J, whereas the
left-hand side tends to -infinity.

Contradiction.

Hence:

    boxed:
    fixed edge + triangle/holonomy couplings do not produce
    p~1/n.

One still needs

    mu_eff ~ log n

or another structural constraint.

## 4. Why holonomy does not rescue sparsity

In the bounded-degree regime, a homogeneous random graph contains too few
active triangles for a fixed triangle-flux term to cancel the combinatorial
edge entropy.

A strong attractive triangle term can instead favor edge condensation into
dense clusters, which is a different phenomenon and does not give homogeneous
finite-degree locality.

## 5. Result

One common gauge-equivariant action can support nontrivial edge phase and
holonomy AFTER active bonds exist.

But in its minimal fixed-coupling form it fails the sparsity gate.

The obstruction is combinatorial, not phase-theoretic.
