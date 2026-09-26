# DYNAMIC-RELATION-LAW-158
## Minimal exchangeable binary bonds are either dense or require system-size-dependent tuning

Date: 2026-09-26

Status:
exact theorem for the minimal permutation-equivariant pair-state Gibbs model.

## 1. Minimal relational degree of freedom

Introduce a binary bond for every unordered pair of candidate units:

    b_xy in {0,1}.

No coordinates or pre-existing graph are supplied.

Take the smallest permutation-equivariant additive relation energy

    E_rel(B)
      =
      epsilon sum_(x<y) b_xy.

The Gibbs/MaxEnt law factorizes over pairs.

Each edge is present with probability

    boxed:
    p
      =
      1/(1+exp(epsilon)).

## 2. Fixed microscopic epsilon gives a dense graph

For n units,

    E[degree]
      =
      (n-1)p.

If epsilon is any fixed finite constant, then p is a fixed nonzero number and

    boxed:
    E[degree] = Theta(n).

Thus the graph becomes dense as n grows.

The only fixed-epsilon sparse limit is

    epsilon -> +infinity,

which gives the empty graph.

## 3. Finite-degree locality requires log(n) tuning

To keep a target mean degree d independent of n, one must impose

    p_n=d/(n-1).

Therefore

    boxed:
    epsilon_n
      =
      log[(n-1-d)/d]
      =
      log(n/d)+o(1).

So the relational chemical potential must know the total system size.

This violates the desired criterion that one fixed microscopic law should work
without a separately tuned graph parameter at every n.

## 4. State-dependent pair energies do not cure the scaling

If each unit has one of q intrinsic types and

    epsilon -> epsilon_ab

for a type pair (a,b), then

    p_ab
      =
      1/(1+exp(epsilon_ab)).

For fixed type fractions x_b, a type-a unit has expected degree

    E[d_a]
      approximately
      n sum_b x_b p_ab.

Again this is O(n) whenever any relevant p_ab is fixed and nonzero.

Thus finite-state typing produces dense block graphs, consistent with report
154.

## 5. Consequence

A minimal binary relational layer DOES create pair-specific random variables,
but fixed local Gibbs weights do not create sparse finite-degree locality.

To get a sparse graph one needs at least one additional mechanism:
- n-dependent chemical potential;
- a hard/soft valence constraint;
- a conservation law for total edge count;
- higher-order relational interactions producing a new scaling regime.

Each is new structure that must itself be sourced.

So DYNAMIC-RELATION-LAW-158 fails in its minimal independent-edge form.
