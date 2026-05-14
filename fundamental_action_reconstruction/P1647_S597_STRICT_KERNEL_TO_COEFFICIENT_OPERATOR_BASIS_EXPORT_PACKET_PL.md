# P1647 / S597 — Strict kernel->coefficient operator-basis export

## Cel
Wzmocnić początek toru (`K_strict -> współczynniki`) przez jawny eksport
mapy projekcyjnej z kernela strict na współczynniki pełnego lagranżianu.

## Zakres
- definicja `K_strict`,
- operatorowa mapa projektorów `I*` do współczynników Lagrangianu,
- jawne zobowiązanie odwracalności (`coefficients -> K_strict`) bez fałszywego PASS.

## Wyjście
- `generated/p1647_s597_strict_kernel_to_coefficient_operator_basis_export_summary.json`
