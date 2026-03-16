# P473 Current Strict Extension‑Lane Global Oriented Selector State — Projector Consistency Audit Probe (No False‑PASS)

Status: `P473_EXECUTED_CURRENT_STRICT_EXTENSION_LANE_GLOBAL_ORIENTED_SELECTOR_STATE_PROJECTOR_CONSISTENCY_AUDIT_PROBE_NO_FALSE_PASS`  
As of: `2026-03-16`

## Goal

After `F470/N516`, strict core exports a **global projective selector state object** on `C_v1` (ray/projector semantics).

After `AX29`, the repo additionally exports an **extension‑lane oriented vector representative** of the same global state by fixing
the residual sign gauge using the directed observable `S_dir(u_1)>0` (`AX28`), applied to the already exported lane‑scoped glued vector section (`F467`).

This probe audits the one no‑false‑pass hygiene condition:

```text
the extension‑lane sign fix changes only a representative vector section,
but does not change the underlying strict projector/ray-level state.
```

## What it checks (declared scope)

For each chart `pair_m (m=1..5)`:

1. Load the extension‑lane oriented vector `u_m` (`AX29` output).
2. Express `u_m` in the declared Fourier chart basis `(c_m,s_m)` and build the 2×2 projector `P_ext = |u><u|`.
3. Load the strict projector operator `A_m(pair_m)` used by the strict projective selector state (`F470` via the exported projector section reference).
4. Check `P_ext` equals `A_m(pair_m)` within a fixed numeric tolerance.

## Exported artifacts

Executed by:

```text
python3 fundamental_action_reconstruction/p473_current_strict_extension_lane_global_oriented_selector_state_projector_consistency_audit_probe.py
```

Exports:

- `fundamental_action_reconstruction/generated/p473_current_strict_extension_lane_global_oriented_selector_state_projector_consistency_audit_probe.json`
- `fundamental_action_reconstruction/generated/p473_current_strict_extension_lane_global_oriented_selector_state_projector_consistency_audit_probe_summary.json`

## Hard limits (no false pass)

`P473` does **not** claim:

1. strict-core discharge of `H37`,
2. strict-core canonical fixing datum export (`T164`),
3. strict-core selector closure / admissible `S_sel_int`,
4. global discharge of `QW-2191`,
5. ToE closure.
