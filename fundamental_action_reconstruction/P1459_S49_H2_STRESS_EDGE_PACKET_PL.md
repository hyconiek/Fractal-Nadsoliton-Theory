# P1459 — S4.9 H2 Stress-Edge Local Replay (PL)

Status: `P1459_EXECUTED_LOCAL_ONLY_STRESS_EDGE_NO_GLOBAL_CLAIM`
As of: `2026-05-13`

## Cel

Wyznaczyć lokalną granicę stabilności remediacji `h2` przez zejście delty w dół
od `base_delta` aż do pierwszego FAIL.

## Procedura

1. Weź baseline z P1452 i `base_delta` z P1453.
2. Schodź deltą: `base_delta, base_delta-step, ...`.
3. Sprawdzaj kryteria:
   - `gain >= min_gain`,
   - `replay_gap <= replay_tol`.
4. Zatrzymaj się na pierwszym FAIL i wyeksportuj obstruction.

## Dyscyplina

- `scope = LOCAL_ONLY_NON_GLOBAL_CLAIM`,
- `strict_core_qw2191_closed = false`,
- `legacy_bridge_used = false`.
