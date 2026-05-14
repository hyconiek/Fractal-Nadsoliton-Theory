# P1631 / S581 — Cover-wise Jacobian invertibility check (strict bidirectional chain)

## Cel
Zrealizować kolejny krok po P1630: sprawdzić na coverze `C_global_noncyclic_cover`
lokalną odwracalność mapy `kernel_params -> coeff` jako warunek dla dwukierunkowego toru
`K_strict <-> coeff <-> L_total <-> EOM`.

## Wejścia
- `generated/p1630_s580_bidirectional_strict_chain_check_summary.json`
- `generated/p1629_s579_global_noncyclic_cover_export_summary.json`

## Wyjście
- `generated/p1631_s581_cover_wise_jacobian_invertibility_summary.json`

## Rygor
- strict-only, bez legacy bridge.
- pełny Lagrangian/EOM jako semantyka toru.
- status global theorem pozostaje OPEN do formalnego dowodu.
