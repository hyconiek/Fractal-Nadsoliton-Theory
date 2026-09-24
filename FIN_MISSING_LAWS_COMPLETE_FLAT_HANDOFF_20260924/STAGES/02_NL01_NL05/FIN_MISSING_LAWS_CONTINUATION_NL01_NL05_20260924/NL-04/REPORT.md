# NL-04 — can refinement/composition/conservation source reversible kinetics?

## Result

`NO_GO_UNIQUE_REVERSIBLE_KINETICS_FROM_CURRENT_AXIOMS`.

The intrinsic spatial pair `(M,K)` from NL-02 admits both

`M qdot = -K q`  (gradient/heat)

and the formal reversible lift

`M qddot + K q = 0`.

Both preserve locality, independent direct-sum composition and the same spatial
refinement laws.  Static/refinement data therefore do not select the dynamical
category.  This recovers the earlier category obstruction by a different route.

Adding energy conservation removes the heat member but still does not make the
reversible theory unique.  NL-03 supplies a four-dimensional internal
deformation fiber.  A kinetic energy

`T = 1/2 sum_i v_i dot(r_i)^T G dot(r_i)`

is compatible with locality, spatial refinement, composition and time reversal
for **any** positive 4x4 tensor G that is used consistently across cells.
Non-proportional G matrices change relative transverse frequencies and cannot be
absorbed into one global clock rescaling.  Hence a new internal kinetic metric
is required.

Even after a kinetic metric is chosen, the current repo still does not source
an independent momentum variable or nondegenerate symplectic form: P3077 found
only formal Hamiltonian lifts, and P3078 accepted zero intrinsic momentum /
symplectic sources.

A further candidate exists inside FIN's statistical structure — Fisher geometry
— but the repo's own ST7665/ST7755/ST7765 conclusion is that information-state
geometry still needs a kinetic metric/reference process to determine motion.
So Fisher geometry does not close NL-04 by itself.

## Minimal missing datum

At minimum one must source:

1. a kinetic metric on the four-dimensional transverse fiber (and on alpha if it is dynamical),
2. an independent momentum/flow or equivalent reversible state extension,
3. a bracket/evolution principle,
4. one overall clock normalization if physical time is intended.
