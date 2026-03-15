# P471 Current Strict `pair1..pair5` Oriented Transport — Operator-Level Cocycle Audit Probe (No False‑PASS)

Status: `P471_DRAFT_EXPECTED_EXECUTED_BY_p471_SCRIPT_NO_FALSE_PASS`  
As of: `2026-03-16`

## Goal

After `F467/P470/N511`, the repo exports a lane-scoped oriented transport lift (`α mod 2π`) on `{pair1..pair5}` as a tracked gauge/convention layer,
with **full triple cocycle/path-independence audited at vector level** on the exported representative vectors.

This probe adds one strict hygiene clarification needed to avoid a common false-pass misread:

```text
the F467 oriented transport cocycle is not an operator-level equality O_jk O_ij = O_ik on the full carrier;
it holds only on the exported glued vector section (and therefore on transported projectors/rays).
```

So `P471` audits, on the exported artifacts:

1. **operator-level** triple residuals: `|| O_jk O_ij - O_ik ||_max`,
2. **vector-section-level** triple residuals: `|| (O_jk O_ij) u_i - (O_ik) u_i ||_2` (consistency with `P470`).

## Inputs

- `F461`: `O_12` (sigma-int corridor).
- `F467`: oriented transport operators `O_13, O_14, O_15, O_23, O_24, O_25, O_34, O_35, O_45`.
- `F456/F462/F464/F465`: exported representative vectors `u_1..u_5` (via `A_1..A_5`).

## Outputs

Executed by:

```text
python3 fundamental_action_reconstruction/p471_current_strict_pair12345_oriented_transport_operator_level_cocycle_audit_probe.py
```

Exports:

- `fundamental_action_reconstruction/generated/p471_current_strict_pair12345_oriented_transport_operator_level_cocycle_audit_probe.json`
- `fundamental_action_reconstruction/generated/p471_current_strict_pair12345_oriented_transport_operator_level_cocycle_audit_probe_summary.json`

## Hard limits

This probe does **not** claim:

- a global selector transition object or global selector atlas on `C_v1`,
- discharge of `QW-2191`,
- any sign-sensitive physical orientation datum (this is a tracked convention layer),
- ToE closure.

