# P1204 W4 Provider Fallback Matrix Packet

Status: `P1204_EXECUTED_W4_PROVIDER_FALLBACK_MATRIX_NO_FALSE_PASS`
As of: `2026-05-10`

## Goal

Następny uczciwy krok po `P1203`: dodać deterministyczny wybór providera CAS
z fallbackiem, aby statusy infrastrukturalne były jednoznaczne.

## Professor-level decision

Dodaję `p1204_w4_provider_fallback_matrix.py`, który:

1. sprawdza dostępność providerów (`sympy`, `sage`, external service),
2. wybiera providera wg priorytetu,
3. aktywuje fallback-mode jeśli żaden provider nie jest dostępny.

## Honest boundary

`P1204` nie wykonuje jeszcze dowodu W4; porządkuje routing backendu.
