# P1456 — S4.6 Semantic Artifact Diff Checkpoint (PL)

Status: `P1456_EXECUTED_LOCAL_ONLY_SEMANTIC_DIFF_NO_GLOBAL_CLAIM`
As of: `2026-05-13`

## Cel

Po technicznym diffie plików (`P1455`) wykonać semantyczny diff local-only dla
małej, krytycznej paczki strict-side (`P1452/P1453/P1455`) bez legacy bridge.

## Zakres semantyczny

- `status` (werdykt),
- `scope_of_pass` lub `scope`,
- `strict_core_qw2191_closed`,
- `legacy_bridge_used`.

## Zasada

Jeśli którykolwiek artefakt utraci `LOCAL_ONLY_NON_GLOBAL_CLAIM` albo ustawi
`strict_core_qw2191_closed=true` bez nowego ścisłego źródła selektora, checkpoint
ma wyjść jako `FAIL_SEMANTIC_GUARDRAIL` i wyeksportować obstruction.

## Dyscyplina

- strict-only,
- non-bridge,
- brak global closure claims.
