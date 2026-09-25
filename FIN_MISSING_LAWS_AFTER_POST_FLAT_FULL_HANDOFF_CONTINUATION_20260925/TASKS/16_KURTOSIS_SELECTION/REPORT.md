# MEMORY-KURTOSIS-07/08 — Fourier selection and exact finite-N witness

Status: **RECOMPUTED_NUMERIC_PLUS_SELECTION_RULE**

## Main result
In uniform equilibrium, the hidden modes are k=1,2 and retained modes are k=3,4,5,6.  Pure quartic self-coupling reaches the hidden sector only for k=5 because `2*5=10=-2 mod 12`.  Exact finite-N generator recursion shows an `O(1/N)` three-generator closure defect for k=5 while pure k=3,4,6 controls are one order smaller / vanish at that leading scale.

## Core formulas
For `f_k=(sqrt(N) mu_k)^4`, define
\[
\Delta_3=(L_{full}^3-L_{ME7}^3)f_k(p_0).
\]
Then k=5 has `N Delta_3 -> nonzero`, while k=3,4,6 have zero leading Fourier overlap.

## Evidence / reproduction
The closed asymptotic formula is recorded in Task 17.  The exact finite-N recursion was used during the research wave; a lightweight formula replay is included here.

## Caveats
The comparison is against the natural max-entropy 7D closure, not every conceivable fitted 7D stochastic model.

## Next question
Use k=3,4,6 as preregistered null controls in a full finite-N simulation.
