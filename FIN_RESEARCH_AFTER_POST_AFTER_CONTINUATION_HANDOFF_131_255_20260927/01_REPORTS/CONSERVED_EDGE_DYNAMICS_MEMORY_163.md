# CONSERVED-EDGE-DYNAMICS-MEMORY-163
## Edge-count conservation does not make neighborhoods persistent under natural reversible swaps

Date: 2026-09-26

Status:
exact one-edge correlation theorem for a reversible fixed-M edge-swap chain.

## 1. Reversible swap dynamics

Keep exactly M active edges among P possible pairs.

Let every active edge ring at rate rho.

When active edge e rings:
- remove e;
- choose uniformly one of the P-M inactive edges;
- activate it.

This preserves M exactly and has the uniform fixed-M ensemble as stationary
law.

## 2. Exact marginal generator

For any fixed candidate edge a:

### if a is active

it disappears at rate

    rho.

### if a is inactive

any of the M active edges may be swapped into a, giving total activation rate

    M rho/(P-M).

Therefore the centered indicator

    b_a-M/P

is an exact eigenfunction.

Its decay rate is

    boxed:
    gamma_edge
      =
      rho P/(P-M)
      =
      rho/(1-p),

where

    p=M/P.

Hence

    boxed:
    Corr[b_a(t),b_a(0)]
      =
      exp[-rho t/(1-p)].

## 3. Sparse limit

For M=dn/2,

    p=d/(n-1)
      ->0.

Therefore

    gamma_edge -> rho

and

    tau_edge -> 1/rho.

So exact edge-count conservation does NOT make individual neighborhoods
long-lived.

At the natural common microscopic clock rho=1, an edge forgets itself on O(1)
time even as n grows.

## 4. Consequence

A conserved global edge count gives:
- sparse instantaneous graphs;

but under ergodic reversible swap dynamics:
- rapidly changing neighbors;
- annealed long-time coupling.

Making neighborhoods persistent requires:
- rho << 1; or
- an additional barrier/local constraint suppressing rewiring.

Choosing rho<<1 by hand reintroduces a free kinetic ratio.

Therefore charge conservation alone does not solve both sparsity and
persistence.
