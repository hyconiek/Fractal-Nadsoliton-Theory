# LOCAL-RELATIONAL-CHARGE-165
## Exact per-unit valence conservation still does not select a unique incidence geometry

Date: 2026-09-26

Status:
exact finite-graph counterexample.

Assume every unit carries an exact local relational charge

    degree(x)=d.

This is much stronger than conserving only the total edge count.

It guarantees finite valence and therefore true graph sparsity when d remains
fixed as n grows.

The question is whether this local conservation law determines WHICH graph is
physical.

## Counterexample: n=6, d=3

Consider two connected 3-regular graphs on six identical units.

### Graph A: K_(3,3)

Every vertex has degree 3.

Laplacian spectrum:

    {0,3,3,3,3,6}.

Diameter:

    2.

### Graph B: triangular prism

Every vertex also has degree 3.

Laplacian spectrum:

    {0,2,3,3,5,5}.

Diameter:

    2.

Both graphs satisfy:
- same n;
- same exact degree at every node;
- connectedness;
- strong vertex symmetry;
- no external coordinates.

Yet they have different collective spectra and different resistance geometry.

Therefore:

    boxed:
    exact valence conservation does not determine incidence.

## Consequence

A local charge d_x solves the number-of-neighbors problem but not the
identity-of-neighbors problem.

Different admissible graphs lead, through the conditional law

    g_r=g-kappa ell_r,

to different collective FIN gains.

So even a derived valence law is insufficient for unique physical predictions.

A further relational invariant must distinguish graph topology.
