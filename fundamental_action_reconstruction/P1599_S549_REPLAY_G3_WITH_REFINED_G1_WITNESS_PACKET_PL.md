# P1599 / S549 — Replay G3 with refined G1 witness

## Cel
Po kroku P1598 wykonać powtórkę finalnej bramki G3, aby sprawdzić,
czy pełny tor strict (`kernel -> współczynniki -> Lagrangian -> EOM`) jest gotowy do domknięcia.

## Wejścia
- `generated/p1596_s546_selector_uniqueness_bridge_theorem_object_summary.json`
- `generated/p1598_s548_strict_g1_witness_refinement_summary.json`
- `generated/p1563_s513_strict_kernel_to_euler_lagrange_eom_export_summary.json`

## Wyjście
- `generated/p1599_s549_replay_g3_with_refined_g1_witness_summary.json`

## Rygor
- Strict-only, bez legacy bridge.
- Bez walidacji zewnętrznej.
- Jawna lista blockerów (`missing_exports`, `missing_witnesses`, `missing_theorems`) przy statusie OPEN.
