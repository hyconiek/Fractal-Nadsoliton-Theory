# P1465 — S4.15 Policy Gate Core Refactor (PL)

Status: `P1465_EXECUTED_LOCAL_ONLY_POLICY_CORE_REFACTOR_NO_GLOBAL_CLAIM`
As of: `2026-05-13`

## Cel

Usunąć duplikację logiki gate przez wspólny rdzeń `p1465_policy_gate_core.py`
wykorzystany przez `P1463` i `P1464`.

## Zakres

- wspólna funkcja decyzyjna `gate_decision(delta, band)`,
- wspólny obiekt pasma `PolicyBand(delta_min, delta_max)`,
- walidacja, że P1463/P1464 dalej dają zgodne wyniki local-only.
