# P1458 — S4.8 H2 Robustness Micro-Sweep (PL)

Status: `P1458_EXECUTED_LOCAL_ONLY_ROBUSTNESS_SWEEP_NO_GLOBAL_CLAIM`
As of: `2026-05-13`

## Cel

Po lokalnym sukcesie P1457 sprawdzić mikro-odporność remediacji `h2` przez
warianty `delta_margin_boost_h2 = base ± ε` bez claimów globalnych.

## Procedura

1. Weź baseline `h1,h2,h3` z P1452.
2. Weź `base_delta` z P1453.
3. Przetestuj mini-siatkę: `{base-0.0005, base, base+0.0005}`.
4. Dla każdego wariantu sprawdź:
   - `gain >= min_gain`,
   - `replay_gap <= replay_tol`.
5. Przy pierwszym FAIL wyeksportuj obstruction.

## Dyscyplina

- `scope = LOCAL_ONLY_NON_GLOBAL_CLAIM`,
- `strict_core_qw2191_closed = false`,
- `legacy_bridge_used = false`.
