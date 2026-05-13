# P1461 — S4.11 Band-Edge Verification Replay (PL)

Status: `P1461_EXECUTED_LOCAL_ONLY_EDGE_VERIFICATION_NO_GLOBAL_CLAIM`
As of: `2026-05-13`

## Cel

Zweryfikować ostrość granicy safety-band z P1460 przez test punktów:

- na krawędziach: `delta=0.0010`, `delta=0.0020`,
- tuż poza: `delta=0.0009`, `delta=0.0021`.

## Kryteria

- `gain >= min_gain`,
- `replay_gap <= replay_tol`,
- obstruction export przy pierwszym FAIL.

## Dyscyplina

- `scope = LOCAL_ONLY_NON_GLOBAL_CLAIM`,
- `strict_core_qw2191_closed = false`,
- `legacy_bridge_used = false`.
