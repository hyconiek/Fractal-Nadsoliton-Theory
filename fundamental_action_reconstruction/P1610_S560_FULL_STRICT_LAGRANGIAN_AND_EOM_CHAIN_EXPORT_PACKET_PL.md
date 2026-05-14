# P1610 / S560 — Full strict Lagrangian and EOM chain export

## Cel
Wyeksportować pełny (nieszkieletowy) strict-only Langrażian oraz jawny łańcuch:
`K_strict -> współczynniki -> L_SM + L_GR + L_mix -> równania ruchu`.

## Wejścia
- `generated/p1558_s508_strict_sm_closed_form_coupling_bundle_and_full_lagrangian_skeleton_summary.json`
- `generated/p1559_s509_gr_strict_curvature_transport_bundle_and_full_lagrangian_update_summary.json`
- `generated/p1562_s512_strict_kernel_to_lagrangian_coefficient_derivation_summary.json`
- `generated/p1563_s513_strict_kernel_to_euler_lagrange_eom_export_summary.json`
- `generated/p1609_s559_strict_kernel_to_full_lagrangian_closure_audit_summary.json`

## Wyjście
- `generated/p1610_s560_full_strict_lagrangian_and_eom_chain_export_summary.json`

## Rygor
- Strict-only, bez legacy bridge i bez transferu ról z kernela legacy.
- Finalny status PASS tylko gdy closure audit jest CLOSED i brak missing items.
- Bez walidacji zewnętrznej.
