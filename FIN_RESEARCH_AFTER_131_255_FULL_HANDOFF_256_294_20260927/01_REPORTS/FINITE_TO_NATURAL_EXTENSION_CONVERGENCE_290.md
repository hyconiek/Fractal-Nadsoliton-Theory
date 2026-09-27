# FINITE-TO-NATURAL-EXTENSION-CONVERGENCE-290
## Finite cyclic Markov bridges converge locally to the natural extension over the same base process

Date: 2026-09-27

Status:
- exact transfer-matrix theorem for any primitive finite Markov kernel;
- exact equality for the full-reset kernel U on every non-wrapping finite window;
- this is process convergence over the base, not uniqueness of internal carrier geometry.

Let P be a primitive stochastic matrix on a finite state set S, with stationary law pi.

The stationary two-sided Markov natural extension has cylinder probabilities

    mu_inf(x_0,...,x_r)
      =
    pi(x_0)
    product_(j=0)^(r-1)
    P(x_j,x_(j+1)).

## 1. Finite cyclic bridge

For cycle length L define a probability measure on

    (x_0,...,x_(L-1))

by

    boxed:
    mu_L(x_0,...,x_(L-1))
      =
    [1 / tr(P^L)]
    product_(j=0)^(L-1)
    P(x_j,x_(j+1)),

with

    x_L=x_0.

Cyclic shift of the L records is an invertible deterministic transformation preserving mu_L.

So this is a finite reversible carrier.

## 2. Exact finite-window marginal

For a block of r+1 consecutive records with r<L,

    boxed:
    mu_L(x_0,...,x_r)
      =
    [ product_(j<r) P(x_j,x_(j+1)) ]
    [ (P^(L-r))_(x_r,x_0) ]
    /
    tr(P^L).

For primitive P,

    P^m
      ->
    1 pi

and

    tr(P^L)
      ->
    1.

Therefore for every fixed r,

    boxed:
    mu_L cylinder
      ->
    mu_inf cylinder

as

    L -> infinity.

This is a concrete finite-to-natural-extension convergence criterion over the declared base process.

## 3. Full-reset Q3 event kernel

For the full reset event kernel

    U_ij=1/3,

we have

    U^m=U

for every m>=1 and

    tr(U^L)=1.

Hence whenever

    L>r,

    boxed:
    mu_L(x_0,...,x_r)
      =
    (1/3)^(r+1)

EXACTLY.

So a finite cyclic register of iid uniform trit records reproduces every finite event-time cylinder of the Bernoulli natural extension exactly until the observation window wraps around the register.

## 4. Generic numerical replay

For a non-rank-one primitive reversible 3-state fixture, total-variation distance of a 4-record block (r=3) was:

    L=4:
      0.316408

    L=5:
      0.174069

    L=8:
      0.031266

    L=12:
      0.003925

    L=24:
      8.49e-6

    L=48:
      4.02e-11.

This is the expected exponential local convergence.

## 5. What this does NOT select

Cylinder convergence only tests the base path process.

Two finite carriers related by an internal measure-preserving recoding can have the same projected cylinder statistics while different internal factor interpretations are used.

Therefore:

    boxed:
    convergence to the natural extension over the base
    is necessary for a faithful finite reversible approximation,
    but does not by itself source spatial geometry.

Operational access to the internal record factors remains an additional structure.

## Verdict

P0-290 PASSES.

There is now an exact process-level sense in which the finite cyclic-record construction approximates the canonical natural extension.

For Q3 full reset the approximation is locally exact before wraparound.
