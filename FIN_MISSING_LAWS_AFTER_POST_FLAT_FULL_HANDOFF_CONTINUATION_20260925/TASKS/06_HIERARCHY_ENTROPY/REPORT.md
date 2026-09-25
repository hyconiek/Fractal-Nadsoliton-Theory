# HIERARCHY-ENTROPY-03 — exact colored-tree entropy rate

Status: **PROOF_GRADE_RECURRENCE_PLUS_RECOMPUTED_NUMERIC**

## Main result
For all K-colored perfect binary hierarchies with unordered children, `A_0=K` and `A_(L+1)=A_L(A_L+1)/2`.  The entropy density `h_K=lim 2^-L log A_L` exists.  Recomputed asymptotic values are h4≈0.8308394326, h8≈1.4522234715, h12≈1.8350061780.

## Core formulas
\[
A_{L+1}=rac{A_L(A_L+1)}2,
\]
\[
h_K=\log(K/2)+\sum_{\ell\ge0}rac{\log(1+A_\ell^{-1})}{2^{\ell+1}}.
\]

## Evidence / reproduction
`REPLAYS/replay_hierarchy_entropy.py` regenerates the values from the recurrence.

## Caveats
Exact entropy for a fixed composition vector is a separate count; the all-color rate is an upper envelope.  Equal-proportion asymptotics require composition-control arguments.

## Next question
Use composition-resolved partition recursions in thermodynamic scaling.
