# P1168 Asymmetric Term Semantic Admissibility Packet

Status: `P1168_EXECUTED_ASYMMETRIC_TERM_SEMANTIC_ADMISSIBILITY_NO_FALSE_PASS`
As of: `2026-05-10`

## Goal

Po `P1167` wykonujemy krok semantyczny: sprawdzić, czy zwycięski kandydat
asymetrycznego członu (`sigma,kappa`) jest dopuszczalny jako strict-side
premise, a nie tylko numeryczny fit.

## Professor-level decision

Wprowadzam checklistę semantyczną dla
`A(d)=sigma*(1-exp(-kappa d))`:

1. `sigma >= 0`,
2. `kappa > 0`,
3. `A(0)=0`,
4. `0 <= A(d) <= sigma` (bounded),
5. monotoniczność `A(d)` dla `d>=0`.

## Result

Dla top-kandydata z `P1166`:

- `semantic_admissible = true`.

Czyli kandydat przechodzi warstwę semantycznej dopuszczalności strict-side.

## Artifacts

- script:
  `p1168_asymmetric_term_semantic_admissibility.py`
- summary:
  `generated/p1168_asymmetric_term_semantic_admissibility_summary.json`

## Honest boundary

`P1168` nie jest dowodem domknięcia teorii i nie rozwiązuje `QW-2191`.
To wyłącznie walidacja semantyczna klasy kandydata.
