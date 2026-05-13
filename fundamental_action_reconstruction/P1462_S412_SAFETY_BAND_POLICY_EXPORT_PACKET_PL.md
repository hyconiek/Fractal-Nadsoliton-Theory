# P1462 — S4.12 Safety-Band Policy Export (PL)

Status: `P1462_EXECUTED_LOCAL_ONLY_POLICY_EXPORT_NO_GLOBAL_CLAIM`
As of: `2026-05-13`

## Cel

Utrwalić jawny policy-band dla `h2` i dodać walidator, który blokuje reruny
wychodzące poza pasmo bez nowego kroku aktualizacji polityki.

## Policy

- `delta_min = 0.0010`
- `delta_max = 0.0020`
- `scope = LOCAL_ONLY_NON_GLOBAL_CLAIM`
- `strict_core_qw2191_closed = false`
- `legacy_bridge_used = false`

## Zasada egzekucyjna

Każda propozycja delty poza pasmem musi dawać `FAIL_POLICY_BAND_VIOLATION` i
obstruction export.
