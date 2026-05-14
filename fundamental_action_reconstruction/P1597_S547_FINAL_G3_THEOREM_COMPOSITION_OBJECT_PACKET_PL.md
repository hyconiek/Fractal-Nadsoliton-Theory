# P1597 / S547 — Final strict G3 theorem composition object

## Cel
Zbudować formalny obiekt końcowego twierdzenia G3 dla toru:
`K_strict -> współczynniki -> pełny Lagrangian -> równania ruchu`,
na ścieżce `F_Nadsoliton => L_SM + L_GR`.

## Wejścia
- `generated/p1596_s546_selector_uniqueness_bridge_theorem_object_summary.json`
- `generated/p1588_s538_g1_full_domain_selector_gap_discharge_summary.json`
- `generated/p1563_s513_strict_kernel_to_euler_lagrange_eom_export_summary.json`

## Wyjście
- `generated/p1597_s547_final_g3_theorem_composition_object_summary.json`

## Rygor
- Strict-only (bez legacy bridge).
- Bez walidacji zewnętrznej.
- Domknięcie tylko przy spełnieniu wszystkich bramek (bridge + G1 + EOM).
