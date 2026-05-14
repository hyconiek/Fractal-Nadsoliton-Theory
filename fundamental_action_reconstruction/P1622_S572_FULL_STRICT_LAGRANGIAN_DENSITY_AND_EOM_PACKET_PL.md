# P1622 / S572 — Full strict Lagrangian density and EOM export (no legacy bridge)

## Cel
Wyprowadzić jawny, sektorowy zapis pełnego lagranżianu strict:
`F_Nadsoliton => L_SM + L_GR + L_mix`,
oraz odpowiadające równania ruchu, bez bridge do legacy.

## Wejścia
- `generated/p1562_s512_strict_kernel_to_lagrangian_coefficient_derivation_summary.json`
- `generated/p1563_s513_strict_kernel_to_euler_lagrange_eom_export_summary.json`

## Wyjście
- `generated/p1622_s572_full_strict_lagrangian_density_and_eom_summary.json`

## Rygor
- Strict-only; `legacy_bridge_used=false`.
- Brak deklaracji closure bez theorem-level domknięcia QW-2191.
- Jawne listy brakujących eksportów/witnessów/theoremów.

## Kontrakt naukowy
1. `K_strict(omega,phi,beta,eta)` mapuje na współczynniki efektywne.
2. Współczynniki mapują na pełny `L_total = L_strict + L_SM + L_GR + L_mix`.
3. EL-EOM są eksportowane sektorowo (`psi`, `H`, `A_mu^a`, `g_mu_nu`).
4. Final strict-core closure pozostaje `OPEN`, dopóki nie ma obiektów theorem-level.
