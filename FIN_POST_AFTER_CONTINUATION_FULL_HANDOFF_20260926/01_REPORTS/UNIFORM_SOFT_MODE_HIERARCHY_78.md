# UNIFORM-SOFT-MODE-HIERARCHY-78
## The strict spectrum predicts four distinct symmetry classes of uniform-state bifurcation

Date: 2026-09-26

Status:
- exact Hessian thresholds for the supplied strict spectrum;
- exact lowest-order invariant classification from D12 representation theory;
- no claim that every nonlinear branch has been globally exhausted.

At the uniform state the retained Fourier sectors diagonalize the dual Hessian:

h_k(g)=1/g-lambda_k/12.

Therefore the softening thresholds are:

k=6: g6 = 5.123427551399
k=5: g5 = 5.220554796949
k=4: g4 = 5.455614632676
k=3: g3 = 6.118057519136.

The order is therefore

k6 -> k5 -> k4 -> k3.

These are not equivalent nonlinear events.

### k=6
A one-dimensional sign representation. Odd powers are forbidden. The quartic coefficient is positive:
D4Phi = 0.0761918988037 > 0.
So the primary event is Z2/pitchfork-like.

### k=5
A faithful two-dimensional D12 representation. The first radial nonlinear term is r^4, while angular anisotropy first appears only at degree 12 through r^12 cos(12 phi). The quartic coefficient is positive:
D4Phi = 0.0550374041046 > 0.
Hence the primary instability is nearly O(2)-radial very close to threshold, which explains the extremely small angular Hessian eigenvalues seen numerically near g≈5.22055.

### k=4
A two-dimensional representation factoring through D3. The cubic invariant Re(z^3) is allowed and nonzero, so amplitudes can scale linearly with parameter distance and subcritical folds are natural. Report 77 is the exact scalar realization.

### k=3
A two-dimensional representation factoring through D4. Cubic invariants are forbidden, but quartic anisotropy Re(z^4) is allowed at the same order as the radial quartic.

Main conclusion:

There is no single “uniform FIN bifurcation.” The local emergence law is determined jointly by:
1. which strict Fourier sector softens,
2. the representation carried by that sector,
3. the first symmetry-allowed nonlinear invariant.

Next atom: K5-NEAR-UNIFORM-CENTER-MANIFOLD-79, aimed at deriving the first nonzero angular term and explaining the reflection/generic branch cluster near g≈5.22055.
