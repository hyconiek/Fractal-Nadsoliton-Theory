# P674 Current Strict Global Directed Selector Closure Output Sign Mismatch Audit Probe

Status: `P674_CURRENT_STRICT_GLOBAL_DIRECTED_SELECTOR_CLOSURE_OUTPUT_SIGN_MISMATCH_AUDIT_PROBE_SPEC_NO_FALSE_PASS`  
As of: `2026-03-17`

## Goal

Audit whether the current exported strict **directed** selector state datum on `C_v1` (vector-level lift; premise-based via `T164/T171`)
induces a **chart-independent directed closure outcome** at the level of the exported global output channels `Y_sel(pair_m)`.

This probe is intentionally narrow:

- it does **not** attempt to define a new directed closure object,
- it does **not** claim strict-core selector closure,
- it does **not** claim global kernel-alone `QW-2191` discharge,
- it only checks whether a *directed* (sign-sensitive) **output vector** is globally consistent across charts under the currently exported fixed bases.

## Inputs

- `generated/selector_state_global_c_v1_directed_strict_v1.json` (directed state; premise-based lift),
- `generated/selector_output_operator_global_c_v1_seed_v1_promoted_strict_v1.json` (global promoted output channels `Y_sel(pair_m)`).

## Acceptance (no false pass)

The probe must:

1. compute the directed selector axis coordinates in each pair chart basis `(c_m,s_m)`,
2. push them forward by the exported output channel `Y_sel(pair_m)` into `(o_+,o_-)`,
3. test whether the sign of the `o_+` component is consistent across `{pair1..pair5}` (up to a single global sign flip),
4. if a mismatch exists, report it as a **directed global closure obstruction** under the current export set.

