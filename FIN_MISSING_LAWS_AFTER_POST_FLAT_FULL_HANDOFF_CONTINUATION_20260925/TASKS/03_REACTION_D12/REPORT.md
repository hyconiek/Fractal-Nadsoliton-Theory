# REACTION-COORD-03 — D12-universal scalar transition sector

Status: **PROOF_GRADE_SYMMETRY_PLUS_PREDECESSOR_NUMERICS**

## Main result
D12 acts orthogonally on retained 7D space and preserves the contracted potential and heat-bath mobility.  Therefore a constrained reaction valley for one localized seed is carried to all 12 localized minima by the D12 action; the same scalar functions `U_eff(s)`, `g_s(s)`, and `m_eff(s)` apply to every orbit copy.

## Core formulas
\[
U_g(R_a\mu)=U_g(\mu),\qquad M(R_a\mu)=R_aM(\mu)R_a^T,
\]
\[
\mu_*^{(a)}(s)=R_a\mu_*(s).
\]

## Evidence / reproduction
Orthogonal D12 action follows from `fin_rank7_followup/src/model.py`; numerical 1D accuracy itself remains a predecessor-package record.

## Caveats
Does not prove global uniqueness of the constrained-minimum branch.

## Next question
Validate the 1D valley with interval/continuation methods if it becomes central.
