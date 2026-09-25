# KURTOSIS-CLOSURE-09/10 — closed k=5 asymptotic coefficient

Status: **THEOREM_CANDIDATE_WITH_RECOMPUTED_FORMULA**

## Main result
The large-N leading third-generator closure defect for the pure k=5 quartic reduces to one hidden k=2 overlap.  The resulting closed coefficient is linear in g and depends only on lambda5 among the strict retained eigenvalues.

## Core formulas
\[
\lim_{N	o\infty}N(L_{full}^3-L_{ME7}^3)x_5^4
=rac{\lambda_5^2}{144}(g\lambda_5-6)
=rac{\lambda_5^2}{24}(1-2\gamma_5),
\]
\[
\gamma_5=1-g\lambda_5/12.
\]

## Evidence / reproduction
`REPLAYS/replay_kurtosis_formula.py` rebuilds lambda5 from the strict kernel and evaluates the formula at g_eq.  At coexistence it gives about `0.0934539123`; sign change occurs at `g=6/lambda5`.

## Caveats
The compact algebraic derivation should still be exported line-by-line in a future proof artifact; current status is theorem candidate rather than repository-certified theorem.

## Next question
Produce a symbolic finite-sum derivation and independent exact-rational/trigonometric replay.
