# P1602 / S552 — Apply tail correction and replay strict readiness

## Cel
Zastosować kandydat korekty ogona (P1601) do strict selector-source i sprawdzić,
czy po replay toru `kernel -> współczynniki -> Lagrangian -> EOM` rośnie gotowość G1/G3.

## Wejścia
- `generated/p1581_s531_strict_selector_source_samples.csv`
- `generated/p1601_s551_strict_selector_source_tail_correction_candidate_summary.json`
- `generated/p1596_s546_selector_uniqueness_bridge_theorem_object_summary.json`
- `generated/p1563_s513_strict_kernel_to_euler_lagrange_eom_export_summary.json`

## Wyjście
- `generated/p1602_s552_apply_tail_correction_and_replay_summary.json`

## Rygor
- Strict-only, bez legacy bridge.
- Bez walidacji zewnętrznej.
- Jawne metryki G1 oraz gate review dla G3.
