# INFORMATION-PRESERVING-TRANSFORMATION-183
## Exact information preservation forces a permutation; one orbit forces a cycle

Date: 2026-09-26

Status:
exact finite-state theorem.

This tests the working relational idea

    relation -> transformation -> relation

without assuming physical space.

## 1. Information-preserving deterministic update

Let

    T:V->V

be a deterministic map on a finite state/unit set V.

Require that for EVERY probability distribution X on V,

    H(T(X))=H(X),

where H is Shannon entropy.

Then T must be injective.

Proof:
if x!=y but T(x)=T(y), take X uniform on {x,y}. Then

    H(X)=log 2

but T(X) is deterministic and

    H(T(X))=0,

contradiction.

For finite T:V->V, injective implies bijective.

Conversely every bijection merely permutes probabilities and preserves H.

Therefore:

    boxed:
    deterministic exact information preservation
      iff
    T is a permutation.

## 2. Graph of the transformation

Draw one directed edge

    x -> T(x).

Because T is a permutation:
- every vertex has out-degree 1;
- every vertex has in-degree 1.

Hence the transformation graph is a disjoint union of directed cycles.

The underlying undirected valence is exactly two on cycles of length >2.

## 3. One orbit

If T is transitive/indecomposable:

    for any x,y there exists k with T^k(x)=y,

then the permutation has exactly one cycle.

Therefore

    boxed:
    information-preserving deterministic transformation
      +
    transitivity
      ->
    one oriented cycle C_n.

This derives the two-port topology of report 178 from transformation structure
rather than postulating degree two directly.

## 4. Boundary

Transitivity is still an additional condition unless sourced.

And the transformation graph is, at this stage, a graph of STATE
transformation.

Identifying it with simultaneous physical adjacency is a further typed bridge.

Nevertheless this is the first route in the current campaign where valence two
arises from a clear information principle rather than an arbitrary graph
constraint.
