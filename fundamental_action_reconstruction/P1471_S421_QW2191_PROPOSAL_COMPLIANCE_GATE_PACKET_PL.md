# P1471 — S4.21 QW-2191 Proposal Compliance Gate (PL)

Status: `P1471_EXECUTED_QW2191_COMPLIANCE_GATE_LOCAL_ONLY`
As of: `2026-05-13`

## Cel

Dodać automatyczną bramkę zgodności z anty-wzorcami P1470, aby nie powtarzać
nieskutecznych dróg w kolejnych propozycjach QW-2191.

## Zasada

Każda nowa propozycja musi podać flagi `avoids_AP1..AP4`.
Brak zgodności -> `FAIL_QW2191_COMPLIANCE_GATE` + obstruction.

## Rygor

- bez legacy bridge,
- bez strict-core closure claim,
- local-only evidence.
