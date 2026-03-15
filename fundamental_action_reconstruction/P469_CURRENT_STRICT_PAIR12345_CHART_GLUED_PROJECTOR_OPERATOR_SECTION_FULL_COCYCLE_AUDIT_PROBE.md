# P469 Current Strict `pair1..pair5` Chart‑Glued Projector Operator Section — Full Cocycle Audit Probe (No False‑PASS)

Status: `P469_DRAFT_EXPECTED_EXECUTED_BY_p469_SCRIPT_NO_FALSE_PASS`  
As of: `2026-03-15`

## Goal

Audit the upgraded strict lane‑scoped five‑chart projector section/gluing data from `F466`:

- gluing laws (projector‑level):
  - `A_2 = O_12 A_1 O_12^T`,
  - `A_3 = O_23 A_2 O_23^T`, `A_3 = O_13 A_1 O_13^T`,
  - `A_4 = O_34 A_3 O_34^T`, `A_4 = O_24 A_2 O_24^T`, `A_4 = O_14 A_1 O_14^T`,
  - `A_5 = O_45 A_4 O_45^T`, `A_5 = O_35 A_3 O_35^T`, `A_5 = O_25 A_2 O_25^T`, `A_5 = O_15 A_1 O_15^T`,
- cocycle/path‑independence on the **exported projector section** for **all** triple overlaps on `{pair1..pair5}`:
  - for each triple `(i,j,k)`, compare the transported projector along `O_kj O_ji` vs the direct transport `O_ki` (sign‑free).

This is a probe‑level numeric audit on the exported current `n=12` artifacts only.

## Inputs

- `F456`, `F462`, `F464`, `F465`: exported projector operators `A_1..A_5`.
- `F461`, `F464`, `F465`, `F466`: exported chart transport operators including `O_12`, `O_23`, `O_13`, `O_34`, `O_45`, `O_24`, `O_35`, `O_14`, `O_15`, `O_25`.

## Outputs

Executed by:

```text
python3 fundamental_action_reconstruction/p469_current_strict_pair12345_chart_glued_projector_operator_section_full_cocycle_audit_probe.py
```

Exports:

- `fundamental_action_reconstruction/generated/p469_current_strict_pair12345_chart_glued_projector_operator_section_full_cocycle_audit_probe.json`
- `fundamental_action_reconstruction/generated/p469_current_strict_pair12345_chart_glued_projector_operator_section_full_cocycle_audit_probe_summary.json`

## Hard limits

This probe does **not** claim:

- a theorem‑level global selector atlas,
- discharge of `QW-2191`,
- sign‑sensitive physical orientation datum,
- ToE closure.

