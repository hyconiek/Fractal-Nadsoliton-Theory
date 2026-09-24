# POST-01 — FLOW-SOURCE / hydrodynamic rank

Status: **THEOREM_PLUS_NUMERICAL_REPLAY**

For a connected parent graph with incidence matrix `B`, the image of `B` is
the mean-zero vertex space, dimension `q-1`.  Four independent constraints
(dipole x/y and traceless quadrupole cos2/sin2) reduce the allowed divergence
space to dimension `q-5`.

At the edge level:
- cycle space: `E-q+1`;
- moment-allowed edge flows: `E-4`;
- visible quotient: `(E-4)-(E-q+1)=q-5`.

For q=12 this is exactly **7**, independent of the number of parent edges.
On the strict complete graph replay:
- `E = 66`;
- `dim cycle = 55`;
- `dim ker(CB) = 62`;
- quotient = **7**.

The strict Laplacian also has an exact edge-divergence realization
`A = B D B^T` to Frobenius error `6.661e-16`.
Conditioning edge flows on the four angular moments produces `A7` with
Frobenius error `1.751e-14`.

For distinct angular directions, the moment spaces through order n are
trigonometric-polynomial spaces of dimension `min(q,2n+1)`.  Hence q=12 has
the exact filtration 7 -> 5 -> 3 -> 1 -> 0 after moments 2,3,4,5,6.
On regular C12 these layers are k=3,4,5,6.

Scope: this sources a **visible quotient dimension** conditional on the moment
law. It does not select a unique microscopic parent covariance or a dynamics.
