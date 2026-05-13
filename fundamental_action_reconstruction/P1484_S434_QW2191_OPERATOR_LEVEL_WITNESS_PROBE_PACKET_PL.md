# P1484 — S4.34 QW-2191 Operator-Level Witness Probe (PL)

Status: `P1484_EXECUTED_QW2191_OPERATOR_LEVEL_WITNESS_PROBE_LOCAL_ONLY`
As of: `2026-05-13`

## Cel

Wykonać następny uczciwy krok po P1483: przejść z poziomu wag sektorowych
na poziom **jawnych świadków operatorowych** dla toru:

`F(nadsoliton) => L_SM + L_GR` (strict-only, bez legacy bridge).

## Decyzja profesorska

Testujemy minimalny, audytowalny zestaw świadków:

- SM-like witness: dodatnia stabilność amplitudowa w oknie SP1,
- GR-like witness: dodatnia stabilność geometryczna w tym samym oknie,
- MIX witness: człon mieszany musi być mniejszy od obu nośników sektorowych
  i mniejszy od policy cap.

Brak claimu closure; `QW-2191` pozostaje jawnie otwarte.
