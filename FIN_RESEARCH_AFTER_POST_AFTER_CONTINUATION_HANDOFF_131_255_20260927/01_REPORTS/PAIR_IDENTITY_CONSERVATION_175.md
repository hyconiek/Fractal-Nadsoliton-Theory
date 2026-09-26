# PAIR-IDENTITY-CONSERVATION-175
## Exact conservation of adjacency gives persistence only by turning the graph into a superselection sector

Date: 2026-09-26

Status:
exact structural theorem.

Let the full microscopic state be

    (node variables, edge variables B),

where B is the complete adjacency pattern.

Suppose pair identity is exactly conserved:

    boxed:
    L B = 0

for the microscopic generator L.

Equivalently, the state space splits into disconnected sectors

    Omega
      =
      disjoint union_B
      Omega_B,

and dynamics never changes B.

## 1. Persistence solved

Within one sector B:
- every neighbor identity is permanent;
- relation memory is infinite;
- the graph can act as a quenched background for slower node dynamics.

So exact pair-identity conservation solves the persistence problem maximally.

## 2. Source problem becomes sector selection

However the dynamics cannot choose between different B sectors because it
never moves between them.

Therefore the physical graph must be supplied by:
- the initial condition;
- a prior/superselection measure;
- or an earlier symmetry-breaking process outside the conserved dynamics.

The conservation law answers:

    how is a graph preserved?

It does not answer:

    why this graph?

## 3. Enormous degeneracy

With n units and M conserved edges, the number of labeled graph sectors is

    binomial[n(n-1)/2, M].

For sparse

    M=O(n),

this still grows super-exponentially in the sense

    log(number of sectors)
      =
      O(n log n).

Thus exact conservation creates a huge family of dynamically disconnected
geometries.

## 4. Conclusion

Pair-identity conservation is not a source law.

It converts the incidence problem into a graph-sector selection problem.

A successful FIN theory would need either:
- a unique/very small admissible sector selected by additional invariants; or
- a dynamical prehistory that selects B and then freezes it, with the selection
  probabilities themselves derived.
