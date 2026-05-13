# P1457 — S4.7 Mini-holdout A/B Replay (PL)

Status: `P1457_EXECUTED_LOCAL_ONLY_AB_REPLAY_NO_GLOBAL_CLAIM`
As of: `2026-05-13`

## Cel

Wykonać uczciwy mini-holdout replay (`h1,h2,h3`) po remediacji `h2` z P1453,
z obowiązkowym porównaniem A/B:

- A = baseline z P1452,
- B = baseline + punktowy boost dla `h2`.

Bez legacy bridge, bez global closure claims.

## Kryteria

- każdy case musi mieć `gain >= min_gain`,
- każdy case musi mieć `replay_gap <= replay_tol`,
- przy pierwszym FAIL: natychmiastowy obstruction export.

## Dyscyplina

- `scope = LOCAL_ONLY_NON_GLOBAL_CLAIM`,
- `strict_core_qw2191_closed = false`,
- `legacy_bridge_used = false`.
