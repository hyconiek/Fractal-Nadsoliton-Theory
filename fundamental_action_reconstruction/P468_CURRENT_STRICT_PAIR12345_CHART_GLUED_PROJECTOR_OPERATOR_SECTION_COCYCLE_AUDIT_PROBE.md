# P468 Current Strict `pair1..pair5` Chart‑Glued Projector Operator Section — Local Cocycle Audit Probe (No False‑PASS)

Status: `P468_DRAFT_EXPECTED_EXECUTED_BY_p468_SCRIPT_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

Audit the newly exported strict lane‑scoped five‑chart projector section/gluing data from `F465`:

- gluing laws (projector‑level):
  - `A_2 = O_12 A_1 O_12^T`,
  - `A_3 = O_23 A_2 O_23^T`,
  - `A_4 = O_34 A_3 O_34^T`,
  - `A_5 = O_45 A_4 O_45^T`,
  - plus direct projector gluing on long edges: `A_3 = O_13 A_1 O_13^T`, `A_4 = O_24 A_2 O_24^T`, `A_5 = O_35 A_3 O_35^T`,
- local cocycle/path‑independence on the **exported projector section** for adjacent triple overlaps:
  - `O_23 O_12` vs `O_13` (1‑2‑3),
  - `O_34 O_23` vs `O_24` (2‑3‑4),
  - `O_45 O_34` vs `O_35` (3‑4‑5).

This is a probe‑level numeric audit on the exported current `n=12` artifacts only.

## Inputs

- `F456`, `F462`, `F464`, `F465`: exported projector operators `A_1..A_5`.
- `F461`, `F464`, `F465`: exported chart transport operators `O_12`, `O_23`, `O_13`, `O_34`, `O_45`, `O_24`, `O_35`.

## Outputs

Executed by:

```text
python3 fundamental_action_reconstruction/p468_current_strict_pair12345_chart_glued_projector_operator_section_cocycle_audit_probe.py
```

Exports:

- `fundamental_action_reconstruction/generated/p468_current_strict_pair12345_chart_glued_projector_operator_section_cocycle_audit_probe.json`
- `fundamental_action_reconstruction/generated/p468_current_strict_pair12345_chart_glued_projector_operator_section_cocycle_audit_probe_summary.json`

## Hard limits

This probe does **not** claim:

- a theorem‑level global selector atlas,
- discharge of `QW-2191`,
- sign‑sensitive physical orientation datum,
- ToE closure.
