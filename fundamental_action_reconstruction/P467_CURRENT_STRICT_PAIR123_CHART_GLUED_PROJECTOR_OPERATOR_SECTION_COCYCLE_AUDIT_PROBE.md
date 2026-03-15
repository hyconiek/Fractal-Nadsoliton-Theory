# P467 Current Strict `pair1/pair2/pair3` Chart‑Glued Projector Operator Section — Cocycle Audit Probe (No False‑PASS)

Status: `P467_DRAFT_EXPECTED_EXECUTED_BY_p467_SCRIPT_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

Audit the newly exported strict lane‑scoped three‑chart projector section/gluing data from `F464`:

- gluing laws:
  - `A_2 = O_12 A_1 O_12^T`,
  - `A_3 = O_23 A_2 O_23^T`,
  - `A_3 = O_13 A_1 O_13^T`,
- cocycle/path‑independence on the **exported projector section**:
  - `O_23 O_12` transports `A_1` to the same `A_3` as `O_13`.

This is a probe‑level numeric audit on the exported current `n=12` artifacts only.

## Inputs

- `F456` exported `A_1(pair1)` projector operator object.
- `F462` exported `A_2(pair2)` projector operator object.
- `F464` exported:
  - `A_3(pair3)` projector operator object,
  - `O_23`, `O_13` chart transport operators,
  - three‑chart glued section / atlas objects.
- `F461` exported `O_12` chart transport operator object.

## Outputs

Executed by:

```text
python3 fundamental_action_reconstruction/p467_current_strict_pair123_chart_glued_projector_operator_section_cocycle_audit_probe.py
```

Exports:

- `fundamental_action_reconstruction/generated/p467_current_strict_pair123_chart_glued_projector_operator_section_cocycle_audit_probe.json`
- `fundamental_action_reconstruction/generated/p467_current_strict_pair123_chart_glued_projector_operator_section_cocycle_audit_probe_summary.json`

## Hard limits

This probe does **not** claim:

- a theorem‑level global selector atlas,
- discharge of `QW-2191`,
- sign‑sensitive physical orientation datum,
- ToE closure.
