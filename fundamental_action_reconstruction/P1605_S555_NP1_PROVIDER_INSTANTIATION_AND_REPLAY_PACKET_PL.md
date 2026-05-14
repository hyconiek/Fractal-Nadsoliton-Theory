# P1605 / S555 — NP1 provider instantiation and replay

## Cel
Zainstancjować noncyclic provider `NP1` i wykonać replay G1/G3 po stronie strict-only,
aby sprawdzić realny wpływ na domknięcie toru `F_Nadsoliton => L_SM + L_GR`.

## Wejścia
- `generated/p1581_s531_strict_selector_source_samples.csv`
- `generated/p1604_s554_noncyclic_provider_upgrade_for_g1_bridge_summary.json`
- `generated/p1596_s546_selector_uniqueness_bridge_theorem_object_summary.json`
- `generated/p1563_s513_strict_kernel_to_euler_lagrange_eom_export_summary.json`

## Wyjście
- `generated/p1605_s555_np1_provider_instantiation_and_replay_summary.json`

## Rygor
- Strict-only, bez legacy bridge.
- Bez walidacji zewnętrznej.
- Jawny raport G1 + bridge + G3 gate po instancjacji NP1.
