# FIN Local Research — Release 10.47

## Stopped Unlimited Algebra Campaign

**Creator:** Żuchowski, Krzysztof  
**Affiliation:** Independent Researcher — Fractal Information Theory Project  
**ORCID:** 0009-0002-0909-3613  
**Resource type:** Publication — Preprint  
**Version:** 1.0.0  
**Publication date:** 2026-08-02  
**Language:** English  
**License:** CC BY 4.0

This source accompanies the English PDF. The campaign was manually stopped
after P475 ran too long. No calculation remains active.

## Result ledger

| Program | Execution | Scientific status |
|---|---|---|
| P486 | Completed | **[Computer-assisted proof]** Exact orientation, branch, pivot, and nonzero-reference premises on the P473 box |
| P485 | Failed before Gröbner | **[Open]** No causal-axis verdict |
| P487 | Failed before Gröbner | **[Open]** No relation for `L` |
| P475 unlimited | Manually interrupted | **[Open]** No basis and no no-go |

P485 and P487 both attempted to coerce the coefficient `-4*sqrt(2)` into
`QQ`. This is a coefficient-domain implementation error. For
`alpha = sin(pi/8)`, the exact identity

`sqrt(2) = 2 - 4*alpha^2`

shows that the coefficient already belongs to `Q(alpha)`. The repair is
proved but was not executed.

P486 proves `det(C)>0`, a nonzero active `A/B/u` pivot, a nonzero tangent
reference coefficient, the positive branch `125*(L-b)-36*alpha=0`, and
`det(Bmap)=1`. It does not prove that the invariant axis lies in the causal
equal-endpoint plane.

P475 was stopped during its fourteen-variable exact lexicographic Gröbner
calculation after more than twelve hours. It returned neither a basis nor an
elimination polynomial. This is not a proof of algebraic impossibility.

The PDF contains the complete derivation, failure analysis, epistemic ledger,
artifact map, and bounded restart specification. Automatic restart is not
authorized.
