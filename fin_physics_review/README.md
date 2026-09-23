# Post-MP7 planning evidence

The analysis inspected the required repository notes and selected MP7 proofs,
code and certificates. It did not rerun the large campaigns or change AGENTS.md.

`audit.py` checks present package/claim hashes, resolves unpacked predecessor
inputs, checks the submitted DAG and selected global-certificate structure,
and records bounded exact endpoint checks. `evidence_audit.json` preserves the
scope limitations and discrepancies; this is not a full theorem reimplementation.

`validate_outputs.py` checks all required task fields, ID/dependency consistency,
local Markdown links, kill-test IDs, and the unchanged AGENTS hash. It generates
`FIN_PHYSICS_PRIORITY_TASKS.json` from the Markdown source of truth.

From the repository root:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 fin_physics_review/audit.py
PYTHONDONTWRITEBYTECODE=1 python3 fin_physics_review/validate_outputs.py --write
PYTHONDONTWRITEBYTECODE=1 python3 fin_physics_review/validate_outputs.py
```

These commands write only new planning/audit artifacts. If the input registry
or AGENTS.md changes, review the new evidence rather than blindly updating the
expected hash. The research tasks themselves remain NOT_STARTED.

## Master-handoff intake (2026-09-23)

The earlier paragraph describes the post-MP7 planning audit. The later FIN PHY
research handoff has now been assessed in
[MASTER_INTAKE_20260923.md](MASTER_INTAKE_20260923.md). Its scoped accepted
results and provenance limits are reflected in the appended FIN PHY section of
`../AGENTS.md`.
