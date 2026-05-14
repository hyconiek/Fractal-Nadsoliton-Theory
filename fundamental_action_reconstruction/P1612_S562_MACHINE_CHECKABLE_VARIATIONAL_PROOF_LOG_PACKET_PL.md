# P1612 / S562 — Machine-checkable variational proof log

## Cel
Utworzyć maszynowo-sprawdzalny log wyprowadzenia EOM z pełnego strict Langrażianu
na torze `K_strict -> współczynniki -> L_total -> EOM`.

## Wejścia
- `generated/p1611_s561_strict_variational_consistency_theorem_object_summary.json`
- `generated/p1610_s560_full_strict_lagrangian_and_eom_chain_export_summary.json`

## Wyjście
- `generated/p1612_s562_machine_checkable_variational_proof_log_summary.json`

## Rygor
- Strict-only, bez legacy bridge.
- Jawne term-by-term mapowanie składników Lagrangianu do sektorów EOM.
- Bez walidacji zewnętrznej.
