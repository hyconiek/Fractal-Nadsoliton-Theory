# P1566 S516 Full Lagrangian EOM Bundle Export With Parameter Sensitivity Packet (No Legacy Bridge)

Status: `P1566_PROPOSED_FULL_LAGRANGIAN_EOM_SENSITIVITY_PACKET`
As of: `2026-05-14`

## Cel

Wzmocnić fizyczną wiarygodność toru strict:

`kernel strict -> współczynniki -> lagranżian -> EOM`

przez analizę czułości parametrów \(\lambda_{SM}^{(eff)}, \kappa_{GR}^{(eff)}, \epsilon_{mix}^{(eff)}\).

## Decyzja profesorska

Eksportujemy pełny bundle:

1. bazowy \(\mathcal{L}_{total}\) i EOM,
2. warianty \(\pm 5\%\) dla każdego współczynnika,
3. metrykę stabilności znaku i skali residuów.

PASS jeśli struktura EOM i klasa stabilności nie zmienia się jakościowo
w oknie \(\pm 5\%\).

## Omówienie dla laika

To sprawdzenie, czy teoria jest „krucha”.
Jeśli mała zmiana parametrów nie psuje równań, to model jest bardziej wiarygodny.
