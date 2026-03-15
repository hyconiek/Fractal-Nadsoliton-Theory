# P437 Current Strict R14/R15/N477 Canonical Local-Diagonal Sigma Opposite‑Pair Sums Value Evaluation Harness Probe

Status: `P437_EXECUTED_CURRENT_STRICT_R14_R15_N477_CANONICAL_LOCAL_DIAGONAL_SIGMA_OPPOSITE_PAIR_SUMS_VALUE_EVALUATION_HARNESS_PROBE_NO_FALSE_PASS`  
As of: `2026-03-14`

## Goal

`T166/T167/P434` reduce the strict diagonal/local `pair1` `O(2)`‑cut question to six opposite‑pair sums on the
canonical diagonal residual profile:

```text
Sigma_psi0_psi6, Sigma_psi1_psi7, Sigma_psi2_psi8,
Sigma_psi3_psi9, Sigma_psi4_psi10, Sigma_psi5_psi11.
```

Independently, `N477` packages a strict conditional reduction:

```text
under constant-vacuum stationarity (canonical EoM) and vpsi_k ≠ 0,
the canonical diagonal residual entries (and therefore the opposite-pair sums)
admit an explicit Yukawa-free, K_total-numeric form (R14) with explicit m0^2 floor (R15).
```

`P437` is an **evaluation harness** for that reduced formula:

1. it reads candidate numeric inputs for `(vpsi[0..11], g4[0..11], g6[0..11])`,
2. it loads strict-fixed `K_total` (from `R14`) and `m0^2` (from `R15`),
3. it computes:
   - the Yukawa-free diagonal residual profile `d[0..11]`,
   - the six opposite‑pair sums `Sigma_psi{k}_psi{k+6}`,
   - and the induced `F2(d)` summary (for convenience only),
4. it writes a persisted JSON artifact for downstream manual review and/or feeding into `P434`.

This probe makes **no** strict-derived value claim. It does not discharge `T167` and must not be promoted into strict
core unless its input values are themselves strict-derived with explicit provenance.

## Execution

Executed by:

- `fundamental_action_reconstruction/p437_current_strict_r14_r15_n477_canonical_local_diagonal_sigma_opposite_pair_sums_value_evaluation_harness_probe.py`

Input location (intentionally `null` by default):

- `fundamental_action_reconstruction/generated/p437_input_vpsi_g4_g6_candidate.json`

CLI overrides (for reproducible non-default runs without mutating the designated candidate file):

- `--in <path>` (input JSON path),
- `--out-json <path>` (artifact JSON path),
- `--out-summary <path>` (summary JSON path).

Persisted artifacts:

- `fundamental_action_reconstruction/generated/p437_current_strict_r14_r15_n477_canonical_local_diagonal_sigma_opposite_pair_sums_value_evaluation_harness_probe.json`
- `fundamental_action_reconstruction/generated/p437_current_strict_r14_r15_n477_canonical_local_diagonal_sigma_opposite_pair_sums_value_evaluation_harness_probe_summary.json`

## Result (no false pass)

On the current repo state, this harness is expected to be **NOT COMPUTABLE** until numeric values for the input arrays
are supplied (with explicit provenance discipline). `T168` names the missing strict-derived upstream value-provider
target class for exactly those inputs.

When the inputs become available, the computed six sums can be copied into the `P434` designated input object to
evaluate the strict decision target `T166`.

## Hard limits

`P437` must not claim:

1. that any computed values are strict-derived unless their provenance is separately discharged,
2. discharge of `T166/T167` or `QW-2191`,
3. strict-core theta export, strict-core selector closure, or ToE closure.
