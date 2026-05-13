# P1469 — S4.19 QW-2191 SP1 Perturbation Stability (PL)

Status: `P1469_EXECUTED_QW2191_SP1_STABILITY_LOCAL_ONLY`
As of: `2026-05-13`

## Cel

Sprawdzić stabilność lokalnego sygnału A/B z P1468 pod małymi perturbacjami
parametryzacji premisy SP1.

## Kryterium

Dla wszystkich wariantów w mikro-siatce wymagamy `delta_metric_B_minus_A > 0`.
Przy pierwszym odwróceniu znaku: obstruction export.

## Rygor

- bez legacy bridge,
- bez claimu strict-core closure,
- status pozostaje `NON_STRICT_UNLESS_PROVEN_INTERNAL`.
