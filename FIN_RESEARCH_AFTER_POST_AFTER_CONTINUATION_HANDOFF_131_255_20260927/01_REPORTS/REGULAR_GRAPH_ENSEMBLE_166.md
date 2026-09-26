# REGULAR-GRAPH-ENSEMBLE-166
## Maximum entropy at fixed valence produces a sparse ensemble, not a unique local geometry

Date: 2026-09-26

Status:
exact symmetry statement plus finite counterexample.

For identical units with

    degree(x)=d

for every x, maximum entropy over labeled simple graphs gives the uniform
d-regular graph ensemble.

By permutation symmetry, for any fixed pair i!=j,

    boxed:
    P(i~j)=d/(n-1).

Therefore the ensemble-mean adjacency is

    E[A]
      =
      [d/(n-1)](J-I),

and the ensemble-mean Laplacian is

    E[L]
      =
      [d/(n-1)](nI-J).

So after annealing, the theory again looks like an all-to-all mean-field model.

## Quenched versus annealed fork

### Annealed relations

If degree-preserving dynamics explores the regular-graph ensemble, neighbor
identity is forgotten.

The long-time effective interaction is mean-field.

### Quenched relations

If one particular regular graph is frozen, sparse neighborhoods persist.

But then the graph realization is additional initial data.

Different realizations with the same d can have different spectra and
transport.

Thus a fixed-valence law leaves a physical fork:
- annealed but nonlocal in time average;
- quenched but realization-dependent.

Neither is a unique derived geometry.
