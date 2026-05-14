# P1629 / S579 — Global noncyclic cover export for strict selector theorem

## Cel
Zrealizować stage_1 z P1628: wyeksportować jawny globalny niecykliczny cover i overlap graph
jako bazę pod globalny dowód selector uniqueness.

## Wejścia
- `generated/p1628_s578_globalization_program_for_strict_selector_proof_summary.json`

## Wyjście
- `generated/p1629_s579_global_noncyclic_cover_export_summary.json`

## Rygor
- strict-only, bez legacy bridge.
- opiera się na pełnym torze `K_strict -> coeff -> L_total -> EOM`.
- nie tworzy closure; tylko dostarcza formalny obiekt cover.
