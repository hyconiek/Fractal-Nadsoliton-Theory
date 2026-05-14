# P1607 / S557 — Final G3 replay with NP1 witness + bridge upgrade

## Cel
Wykonać finalną bramkę replay G3 z importem `W1606/T1606` i wydać uczciwą decyzję
`closure vs nonclosure` dla strict-core bez legacy bridge.

## Wejścia
- `generated/p1606_s556_np1_witness_and_bridge_upgrade_summary.json`
- `generated/p1563_s513_strict_kernel_to_euler_lagrange_eom_export_summary.json`

## Wyjście
- `generated/p1607_s557_final_g3_replay_with_np1_witness_bridge_summary.json`

## Rygor
- Strict-only.
- Bez walidacji zewnętrznej.
- Jawne kryterium zamknięcia: witness + bridge + EOM.
