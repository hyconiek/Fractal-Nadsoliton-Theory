# P1595 / S545 — Final strict G3 gate (kernel → współczynniki → Lagrangian → EOM)

## Cel
Domknąć (lub uczciwie pozostawić otwarte) strict-core ToE dla toru:

`F_Nadsoliton => L_SM + L_GR`,

przez końcową bramkę G3, opartą wyłącznie o strict-chain i bez bridge do legacy.

## Zakres rygoru
- Brak transferu roli fizycznej z kernela legacy.
- Brak walidacji przez zewnętrzne zespoły.
- Jawna lista braków w trzech klasach: `missing_exports`, `missing_witnesses`, `missing_theorems`.

## Wejścia
- `generated/p1593_s543_focused_first_obligation_discharge_summary.json`
- `generated/p1594_s544_focused_second_obligation_g2_discharge_summary.json`
- `generated/p1562_s512_strict_kernel_to_lagrangian_coefficient_derivation_summary.json`
- `generated/p1563_s513_strict_kernel_to_euler_lagrange_eom_export_summary.json`
- `generated/p1582_s532_strict_selector_uniqueness_theorem_bridge_to_full_lagrangian_summary.json`

## Wyjście
- `generated/p1595_s545_final_g3_attempt_from_g1_g2_summary.json`

## Kryterium domknięcia
PASS tylko gdy jednocześnie:
1. G1 gotowe,
2. G2 gotowe,
3. eksport `kernel -> coefficients` gotowy,
4. eksport `lagrangian -> eom` gotowy,
5. gotowy theorem-bridge `selector uniqueness -> full Lagrangian`.

W przeciwnym razie status pozostaje `KEEP_OPEN_P1595_G3_NOT_READY`.
