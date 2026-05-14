# P1611 / S561 — Strict variational consistency theorem object

## Cel
Zbudować obiekt theorem-level spinający pełny Langrażian strict z równaniami ruchu
w każdym sektorze (`SM`, `GR`, `mix`) na torze:
`K_strict -> współczynniki -> pełny Langrażian -> EOM -> closure`.

## Wejścia
- `generated/p1610_s560_full_strict_lagrangian_and_eom_chain_export_summary.json`
- `generated/p1609_s559_strict_kernel_to_full_lagrangian_closure_audit_summary.json`

## Wyjście
- `generated/p1611_s561_strict_variational_consistency_theorem_object_summary.json`

## Rygor
- Strict-only; bez legacy bridge.
- Obiekt theorem-level musi jawnie raportować brakujące exports/witnesses/theorems.
- Bez walidacji zewnętrznej.
