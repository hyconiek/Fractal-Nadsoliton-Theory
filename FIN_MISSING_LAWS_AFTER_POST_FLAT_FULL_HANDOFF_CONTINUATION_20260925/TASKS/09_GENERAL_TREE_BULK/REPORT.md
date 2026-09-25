# TREE-BULK-02 — arbitrary-tree local realization theorem

Status: **PROOF_GRADE_CONDITIONAL**

## Main result
For an arbitrary rooted hierarchy with positive leaf masses and node stiffnesses increasing down each branch, the hierarchy quadratic form has an exact positive local tree-network realization.  The edge conductance is determined uniquely in the scalar weighted-mean realization class.

## Core formulas
For internal child c of v: \[
g_{vc}=rac{M_c}{\kappa_v^{-1}-\kappa_c^{-1}},\qquad \kappa_c>\kappa_v.
\]
For a leaf child ell: \[
g_{v\ell}=\kappa_vm_\ell.
\]

## Evidence / reproduction
Inductive series/parallel proof; an irregular 5-leaf numerical Schur test gave max error about `5.6e-17`.

## Caveats
Uniqueness is within the declared scalar-edge/weighted-mean realization class.

## Next question
Classify inverse reconstruction and dynamic storage on this general tree.
