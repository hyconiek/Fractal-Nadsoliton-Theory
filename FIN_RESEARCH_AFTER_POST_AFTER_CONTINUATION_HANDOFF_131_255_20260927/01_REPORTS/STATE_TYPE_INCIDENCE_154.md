# STATE-TYPE-INCIDENCE-154
## Unit states alone can generate only block-complete graphs, not generic locality

Date: 2026-09-26

Status:
exact extension of report 153.

Suppose the n units are not all identical but belong to a finite set of q
intrinsic coarse-state types:

    tau_i in {1,...,q}.

Assume the incidence rule depends only on the unit states/types and remains
permutation-equivariant.

Let n_a be the number of units of type a.

The stabilizer group of the type assignment is

    S_(n1) x ... x S_(nq).

This group is transitive on:
- unordered pairs within the same type a;
- ordered type-pairs (a,b) with a!=b.

Therefore the most general symmetric weight matrix has block-constant form

    W_ij
      =
      w_(tau_i,tau_j).

There are at most

    q(q+1)/2

distinct pair weights.

## Consequence

For fixed q as n grows, the graph is a weighted complete block graph.

State typing can create:
- communities;
- complete multipartite structure;
- finite-rank collective modes.

It cannot create a generic sparse geometry in which two units of the same type
have different neighborhoods.

So even the dynamically selected three-state FIN variable does not by itself
supply physical incidence among many simultaneous units.

With q=3 it can at most generate six symmetry-distinct pair weights.

## Required new object

To obtain local neighborhoods one needs a pair-specific relation

    r_xy

that is not reducible to

    (tau_x,tau_y)

alone.

This relation may later become:
- conductance;
- adjacency;
- distance-like data;
- a dynamical bond state.

But it must be an actual degree of freedom or a derived pair observable.

The single-cell state labels are insufficient.
