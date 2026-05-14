# P1613 / S563 — Symbolic Euler-Lagrange identity replay

## Cel
Dodać artefakt symbolic-replay, który jawnie potwierdza tożsamości
Euler-Lagrange dla pełnego strict Langrażianu (SM/GR/MIX).

## Wejścia
- `generated/p1612_s562_machine_checkable_variational_proof_log_summary.json`
- `generated/p1610_s560_full_strict_lagrangian_and_eom_chain_export_summary.json`

## Wyjście
- `generated/p1613_s563_symbolic_euler_lagrange_identity_replay_summary.json`

## Rygor
- Strict-only, bez legacy bridge.
- Tożsamości muszą mieć jawny status `verified` per sektor.
- Bez walidacji zewnętrznej.
