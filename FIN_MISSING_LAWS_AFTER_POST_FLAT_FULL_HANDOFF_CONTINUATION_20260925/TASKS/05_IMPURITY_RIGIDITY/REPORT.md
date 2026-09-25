# DUPLICATE-GAP-02 — impurity-defect scaling and rigidity bound

Status: **PROOF_GRADE_FOR_DECLARED_DEFECT_CLASS**

## Main result
A single impurity propagated through a dyadic replicated subtree has cost growing as `n^(log_2 s)=n^(2/D_H)`.  Therefore a necessary condition for a strictly positive hierarchy-gap density under this defect mechanism is `D_H<=2` (`r>=sqrt(2)` in the binary case).

## Core formulas
For `n=2^m`, \[
E_{m imp}(n)=rac{s^L d^2}{n}rac{(2s)^m-1}{2s-1}
\sim C n^{\log_2 s}=C n^{2/D_H}.
\]


## Evidence / reproduction
Algebraic within the specified impurity construction.

## Caveats
`D_H<=2` is a necessary condition for extensive rigidity in this defect class, not a proof that D_H>2 cannot order by other mechanisms.

## Next question
Compare the full defect spectrum with hierarchy entropy.
