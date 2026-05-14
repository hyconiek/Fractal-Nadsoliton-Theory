# P1638 / S588 — Strict selector-constraint candidate for nullspace removal

## Cel
Po P1637: zaproponować i przetestować kandydat warunku selektora, który usuwa degenerację
`coeff -> kernel` (nullspace) w strict-only torze.

## Założenie robocze
Używamy warunku atlasowego: `eta` jest ustalone przez globalny noncykliczny cover witness
(z P1631/P1632) i nie może dryfować przy inwersji lokalnej.

## Wyjście
- `generated/p1638_s588_strict_selector_constraint_candidate_for_nullspace_removal_summary.json`
