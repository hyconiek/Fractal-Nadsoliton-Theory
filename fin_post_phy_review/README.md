# Post-PHY architecture review

This directory contains the evidence index and the source of the new research
design. It does not execute its research targets or modify AGENTS.md.

- `search_index.json`: repo-wide content-search groups, overlapping file lists,
  and logical source paths/hashes from the master content store.
- `selected_sources.json`: exact payload identities cited by task cards.
- `program_spec.py`: 42 task records with explicit de-dup, methods, claims,
  stops, held-out tests, dependencies and eight-dimensional scores.
- `build_program.py`: renders the catalogue, task JSON, DAG, de-dup table and
  a synchronized extract of kill/source/physical gates.
- `validate_program.py`: validates fields, links, DAG, source identities and
  the unchanged guardrail hash. It is not a mathematical theorem verifier.

The master was hash-checked and materialized to
`/tmp/fin_post_phy_design_20260923` for inspection. No research producer was
rerun as part of this planning task. The science baseline is HEAD `01c89aa3`;
the intervening `d4c4a0ac` commit stores planning work only. Validation records
the actual HEAD and rejects unreviewed intervening scientific changes.

To regenerate the design from the repository root:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 fin_post_phy_review/build_program.py
PYTHONDONTWRITEBYTECODE=1 python3 fin_post_phy_review/validate_program.py
```

`search_audit.py` repeats the broad content search and can be slower on this
large repository. It is not necessary just to render the catalogue. A new
science HEAD requires a new review, not editing the expected hash blindly.

Primary-source checks through Firecrawl informed the use of safe likelihood
upper bounds and circular-planar response criteria. Links and their limited
role are documented in the main programme. The proposed escape–profile,
mesoscopic-window and source/closure bridges remain research targets, not
newly accepted physics.

## Research-handoff intake (2026-09-24)

The planning description above predates execution. The later 42-node handoff is
assessed in [RESEARCH_INTAKE_20260924.md](RESEARCH_INTAKE_20260924.md).
Its accepted scopes and fail-closed provenance limits are now recorded in the
FIN post-PHY section of `../AGENTS.md`.
