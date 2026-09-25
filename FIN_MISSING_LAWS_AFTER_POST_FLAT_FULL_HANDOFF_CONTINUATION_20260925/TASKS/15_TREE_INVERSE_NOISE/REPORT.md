# TREE-INVERSE-NOISE-06 — sufficient topology stability under response noise

Status: **PROOF_GRADE_PERTURBATION_BOUND**

## Main result
If a measured static boundary Laplacian has spectral error epsilon below the smallest positive eigenvalue sigma, the induced effective-resistance error is bounded.  If that error is less than half the shortest internal edge length of the true additive tree, the recovered topology is guaranteed unchanged.

## Core formulas
\[
\delta_R\le rac{2arepsilon}{\sigma(\sigma-arepsilon)}.
\]
A sufficient topology-preservation condition is \[
arepsilon<rac{w_{\min}\sigma^2}{4+w_{\min}\sigma}.
\]

## Evidence / reproduction
Obtained from pseudoinverse perturbation on the mean-zero sector plus the four-point/tree-metric topology margin.

## Caveats
Assumes the zero mode is preserved/aligned and the minimal internal edge resistance `w_min` is positive.

## Next question
Extend to noisy dynamic residues/storages and confidence regions.
