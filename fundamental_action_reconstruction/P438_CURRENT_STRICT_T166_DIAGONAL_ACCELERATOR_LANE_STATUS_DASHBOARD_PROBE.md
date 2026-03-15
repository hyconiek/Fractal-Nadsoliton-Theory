# P438 Current Strict T166 Diagonal Accelerator Lane Status Dashboard Probe

Status: `P438_EXECUTED_CURRENT_STRICT_T166_DIAGONAL_ACCELERATOR_LANE_STATUS_DASHBOARD_PROBE_NO_FALSE_PASS`  
As of: `2026-03-14`

## Goal

The diagonal/local accelerator lane for `pair1` (the `T166` decision target) is now spread across multiple strict
artifacts:

1. repo-wide value-instantiation scan (`P432`),
2. sigma-six-sum evaluator (`P434`),
3. Yukawa-free `N477` sigma-six-sum harness (`P437`),
4. missing-target specs (`T167`, `T168`).

This creates a usability risk:

```text
people keep asking “is it computable yet?” and accidentally re-run or verbally promote toy values.
```

`P438` is an *audit dashboard* probe: it runs the mechanical probes (`P432/P444/P437/P434`) and emits one single JSON
summary capturing:

- what is computable today,
- exactly which numeric inputs are missing,
- and which strict target spec (`T167` vs `T168`) is the next required move.

Additionally, it runs `P444` as a repo-state hygiene check for any already-exported `T168`-consumable `(vpsi,g4,g6)`
value-provider object (without promoting any toy/extension values).

It exports **no new physics** and makes **no** closure claim.

## Execution

Default (fast hygiene scan for `P432`):

```bash
python3 fundamental_action_reconstruction/p438_current_strict_t166_diagonal_accelerator_lane_status_dashboard_probe.py
```

Full repo scan (run `P432` with `--full-scan`, i.e. include vendor/cache directories):

```bash
python3 fundamental_action_reconstruction/p438_current_strict_t166_diagonal_accelerator_lane_status_dashboard_probe.py --full-scan
```

## Inputs reused

1. `P432` — repo scan for decision-ready numeric instantiation at the designated location.
2. `P444` — repo scan for any `T168`-consumable `(vpsi,g4,g6)` instantiation and strict-derived marker presence.
3. `P437` — computes the six opposite-pair sums from `(vpsi,g4,g6)` **if** those inputs are provided.
4. `P434` — computes `F2(d)` and induced `pair1` anisotropy **if** the six sigmas are provided.
5. `T166` — decision target (context).
6. `T167/T168` — missing strict-derived value-provider targets (context).

## Output

`P438` persists:

- `fundamental_action_reconstruction/generated/p438_current_strict_t166_diagonal_accelerator_lane_status_dashboard_probe.json`
- `fundamental_action_reconstruction/generated/p438_current_strict_t166_diagonal_accelerator_lane_status_dashboard_probe_summary.json`

## Hard limits (no false pass)

`P438` does **not** claim:

1. discharge of `T166` (no decision of `F2(d)` for canonical `D_local_residual`),
2. strict-derived value export for the six sigmas (`T167`),
3. strict-derived value export for `(vpsi,g4,g6)` (`T168`),
4. strict-core selector closure or `QW-2191` discharge,
5. ToE closure.
