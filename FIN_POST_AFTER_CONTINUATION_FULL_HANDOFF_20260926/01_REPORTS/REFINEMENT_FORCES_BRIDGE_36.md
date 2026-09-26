# REFINEMENT-FORCES-BRIDGE-36
## Exact edge-refinement rigidity, and the typed barrier to the strict intracell operator

Date: 2026-09-26

Repository baseline:
`hyconiek/Fractal-Nadsoliton-Theory`, visible HEAD
`fe14a6f4e436815635df54429102f23c22862296`.

Status:
- **EXACT theorem** in the already admitted metric-edge Dirichlet refinement lane;
- **EXACT typed no-go** for transferring that theorem directly to the full dense
  strict C12 operator;
- no new physical units, clock or GR claim.

This continues `DYNAMIC-HIDDEN-BRIDGE-SOURCE-35`.

## 1. Accepted refinement input

The existing FIN geometry/source work gives, conditionally on the declared
state-sourced circle geometry and exact static subdivision principle,

    c(ell)=kappa0/ell.

Thus resistances add:

    R(ell)=ell/kappa0.

For an interval split as

    ell=a+b,

static Schur/Kron elimination of the midpoint reproduces exactly c(ell).

## 2. Hidden scalar interpolation

Let a scalar hidden/modulating field have endpoint values u0,u1.

If the same Dirichlet law is used to eliminate the inserted midpoint, its
harmonic value is

    um = (b u0 + a u1)/(a+b).

This is just the linear interpolant on the metric interval.

## 3. General linear symmetric fractional modulation

Consider a first-order fractional conductance modulation on a segment of
length ell:

    c_eps(ell;u,v)
      = c(ell) [1 + eps beta(ell)(u+v)/2] + O(eps^2).

The assumptions are:

1. first order in the hidden scalar;
2. symmetric in the two endpoints;
3. no extra internal coordinate beyond segment length;
4. same law is used before and after subdivision.

No desired beta is inserted.

## 4. Exact series-elimination formula

For the two fine segments,

    c_a = kappa0/a [1+eps phi_a] + O(eps^2),
    c_b = kappa0/b [1+eps phi_b] + O(eps^2),

where

    phi_a = beta(a)(u0+um)/2,
    phi_b = beta(b)(um+u1)/2.

Because resistances add, the effective fractional modulation is

    phi_eff
      = [a phi_a + b phi_b]/(a+b).

Exact refinement compatibility requires

    phi_eff
      = beta(a+b)(u0+u1)/2

for every positive a,b and arbitrary u0,u1.

## 5. Rigidity theorem

Substituting the harmonic midpoint and matching separately the coefficients
of u0 and u1 gives

    beta(a)=beta(a+b),
    beta(b)=beta(a+b).

Hence for arbitrary positive a,b,

    beta(a)=beta(b)=beta(a+b).

Therefore beta is constant on positive lengths:

    beta(ell)=beta0.

So exact arbitrary-split refinement forces

    boxed:
    delta c / c
      = eps beta0 (u_left+u_right)/2.

The endpoint-average *fractional* bridge is therefore unique up to one global
coupling beta0 in this metric-edge lane.

This is stronger than the six-shell freedom found in report 35.

## 6. Operator form on an arbitrary metric graph

For a weighted graph with edge conductances wij, define u on vertices.

The first-order perturbation

    delta w_ij
      = beta0/2 * w_ij (u_i+u_j)

induces the Laplacian variation

    delta A
      = beta0/2 [
          D_u A + A D_u - diag(Au)
        ].

Thus the algebraic bridge from report 35 is exactly the graph-level form of
the refinement-rigid edge law.

Within the metric-edge lane, the previously "candidate" bridge is therefore
derived up to one scalar.

## 7. Why this does NOT yet source the bridge for the full strict A

The repository separately proves:

> exact static elimination of purely local subdivisions of a one-dimensional
> coarse cycle can generate only nearest-neighbour coarse couplings.

But the frozen strict C12 operator has nonzero shells d=2,...,6.

The accepted strict non-nearest exit-weight fraction is about

    0.4338569991.

Therefore the full dense strict A is not the exact static coarse Schur image of
that local metric-circle subdivision class.

Consequently the refinement theorem above cannot be applied to the six strict
intracell shells merely by calling d a physical length.

This would cross the open legacy/strict or intracell/intercell bridge without
a theorem.

## 8. Precise outcome of REFINEMENT-FORCES-BRIDGE-36

Positive result:

    metric edge
    + c(ell)=kappa0/ell
    + harmonic hidden interpolation
    + arbitrary split consistency

forces

    delta c/c = beta0 endpoint_average(u)

with one global beta0.

Negative result:

the existing theorem acts on the state-sourced sparse intercell metric graph,
not on the dense strict intracell C12 operator.

Therefore refinement reduces the constitutive freedom to one scalar **after a
typed metric-edge identification exists**, but does not itself provide that
identification for strict A.

## 9. Physical significance

This narrows the missing physics law substantially.

We no longer need to ask:

    "Why this particular endpoint-average conductance formula?"

Inside the metric/refinement lane it is forced.

The remaining questions are instead:

1. what FIN object supplies the scalar field u on each emergent spatial cell?
2. why/how does that field couple to the intercell Dirichlet energy?
3. what fixes beta0 and its dimensional normalization?
4. can the four-dimensional hidden H4 sector supply u without an external
   selector?

The fourth question is attacked in the companion continuation
`PHASE_LOCKED_HIDDEN_SCALAR_37`.

## 10. Boundary

Not concluded:
- full strict A is a physical spatial Laplacian;
- H4 is a physical field;
- beta0 is known;
- metre/second/action units are derived;
- GR backreaction is derived;
- QW-2191, role-bearing L_total, SM/GR or ToE are closed.
