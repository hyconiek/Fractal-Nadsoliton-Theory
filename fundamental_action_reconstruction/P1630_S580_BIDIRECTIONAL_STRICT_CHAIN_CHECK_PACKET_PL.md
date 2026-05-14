# P1630 / S580 — Bidirectional strict chain check

## Cel
Sprawdzić dwukierunkowość toru strict:
- forward: `K_strict -> coeff -> L_total -> EOM`
- backward (lokalny): `coeff -> (omega,phi,beta,eta)`

bez legacy bridge.

## Wejścia
- `generated/p1562_s512_strict_kernel_to_lagrangian_coefficient_derivation_summary.json`
- `generated/p1629_s579_global_noncyclic_cover_export_summary.json`

## Wyjście
- `generated/p1630_s580_bidirectional_strict_chain_check_summary.json`

## Rygor
- strict-only.
- pełny Lagrangian i EOM użyte jawnie.
- backward check oznaczony jako lokalny/inwersyjny, bez nadpisywania statusu global theorem.
