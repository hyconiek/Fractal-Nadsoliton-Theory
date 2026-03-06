# RAPORT QW-2123: VACUUM BRANCH SELECTION STRICT GATE

- Date UTC: 2026-03-04T18:48:07.994011+00:00
- Verdict: **VACUUM_BRANCH_SELECTION_STRICT_GATE_PASS_BROKEN_BRANCH_REQUIRED**
- pass_count: `10/10`

## Inputs
- lambda_min(K_total): `-0.681874762380`
- required_shift >= `0.681874763380`
- symmetric_floor: `0.506775986179`
- broken_floor: `1.013551972358`

## Selection rule
- If lambda_min(K_total)<0 and symmetric_floor<required_shift<=broken_floor, select broken branch as strict physical closure branch.
- result: `broken_branch_required`

## Artifact
- JSON: `report_qw2123_vacuum_branch_selection_strict_gate.json`
