# P1463 — S4.13 Policy-Aware Rerun Gate (PL)

Status: `P1463_EXECUTED_LOCAL_ONLY_POLICY_GATE_NO_GLOBAL_CLAIM`
As of: `2026-05-13`

## Cel

Zintegrować policy-band P1462 jako bramkę wejściową dla rerunu, tak aby
out-of-band delta była blokowana przed uruchomieniem replay.

## Zasada

- jeśli `delta ∈ [delta_min, delta_max]` → `PASS_GATE_IN_BAND`,
- jeśli `delta` poza pasmem → `FAIL_POLICY_BAND_VIOLATION` + obstruction.

## Dyscyplina

- `scope = LOCAL_ONLY_NON_GLOBAL_CLAIM`,
- `strict_core_qw2191_closed = false`,
- `legacy_bridge_used = false`.
