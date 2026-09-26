# STATIONARY-BRANCH-GRAPH-COMPLETION-76
## Every orbit in the numerical g=5 atlas is now assigned to a connected branch component

Date: 2026-09-26

Repository baseline:
`hyconiek/Fractal-Nadsoliton-Theory`, visible HEAD
`fe14a6f4e436815635df54429102f23c22862296`.

Status:
- branch pairing below is NUMERICAL pseudo-arclength continuation;
- previously certified local nodes remain certified at their own locations;
- this is not a global stationary-exhaustion theorem.

## Starting atlas

R7P-037 found 15 D12 orbits at exact g=5.

Earlier continuation reports already assigned:
- i=0: main localized minimum, index 0
- i=2: secondary saddle, index 1
- i=3: large two-harmonic parent, index 2
- i=7: main localization saddle, index 1
- i=9: second D3 daughter, index 2
- i=10: small two-harmonic/D3 parent, index 3
- i=6: uniform branch.

This left eight nonuniform atlas orbits unassigned.

## Four missing components

Pseudo-arclength continuation pairs the remaining eight orbits as follows.

### Component A
i=1, index 1
<-> fold at g ≈ 4.352214678656
<-> i=8, index 2.

Fold spectrum begins approximately:
(-0.1223270, 0, +0.0951751, ...).

### Component B
i=4, index 1
<-> fold at g ≈ 4.890011279965
<-> i=13, index 2.

Fold spectrum begins approximately:
(-0.1229342, 0, +0.0104816, ...).

### Component C
i=5, index 2
<-> fold at g ≈ 4.913081079101
<-> i=14, index 3.

Fold spectrum begins approximately:
(-0.1008894, -0.0147627, 0, +0.0524537, ...).

### Component D
i=11, index 2
<-> fold at g ≈ 4.993057757778
<-> i=12, index 3.

Fold spectrum begins approximately:
(-0.0454285, -0.0454285, 0, +0.00509624, ...).

The last component lies in the exact pure-k=4 invariant family.

## Complete g=5 assignment

All 15 known R7P-037 orbits now have a branch-role assignment:

i=0 main minimum component
i=1 component A
i=2 secondary D3 component
i=3 two-harmonic parent component
i=4 component B
i=5 component C
i=6 uniform branch
i=7 main saddle component
i=8 component A
i=9 secondary D3 component
i=10 two-harmonic/D3 parent component
i=11 pure-k4 component D
i=12 pure-k4 component D
i=13 component B
i=14 component C.

This does not prove stationary exhaustion. It proves only that every orbit already found at g=5 is now placed on a connected branch component.

## Near-uniform clue

Components B and C move toward the soft k=5 sector near the exact uniform threshold

g5 = 12/lambda5 ≈ 5.220554796949.

The reflection-symmetric branch approaches the primary k=5 regime, while the generic component develops an additional very small angular Hessian eigenvalue nearby, suggesting secondary symmetry breaking rather than an independent primary uniform bifurcation.
