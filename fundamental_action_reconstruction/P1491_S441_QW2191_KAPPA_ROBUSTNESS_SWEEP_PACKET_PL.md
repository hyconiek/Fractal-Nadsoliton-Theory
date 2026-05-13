# P1491 — S4.41 QW-2191 Kappa Robustness Sweep (PL)

Status: `P1491_EXECUTED_QW2191_KAPPA_ROBUSTNESS_SWEEP_LOCAL_ONLY`
As of: `2026-05-13`

## Cel

Wykonać następny uczciwy krok do domknięcia QW-2191:
sprawdzić, czy redukcja luki selektora utrzymuje się w całym dopuszczalnym
oknie `kappa`, a nie tylko punktowo.

## Decyzja profesorska

Dla siatki `kappa` testujemy:

1. `G1(kappa) < G0` (realna redukcja luki),
2. `|Delta_SB(kappa)| <= safety_margin` (fizyczna kontrola),
3. monotoniczny trend poprawy w bezpiecznym zakresie.

Jeśli większość bezpiecznych punktów przechodzi, mamy silniejszą podstawę
fizyczną pod próbę domknięcia theorem-level.
