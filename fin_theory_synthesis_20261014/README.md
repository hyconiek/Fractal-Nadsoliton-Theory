# FIN theory-only continuation after GEN and RG-F

Read [RESEARCH_REORIENTATION_PL.md](RESEARCH_REORIENTATION_PL.md), then
[GEN_EVIDENCE_MAP.md](GEN_EVIDENCE_MAP.md),
[CORE_DERIVED_BRIDGES.md](CORE_DERIVED_BRIDGES.md) and
[LITERATURE_DATABASE.md](LITERATURE_DATABASE.md).

Machine records: LITERATURE_REGISTER.json, CLAIM_REGISTER.json,
verification.json and GEN_HEADINGS_INDEX.json. Index counts are narrative
navigation, not numbers of proved results.

Execution plan:
[FIN_INTRINSIC_THEORY_NEXT_CAMPAIGN_20261014.md](../FIN_INTRINSIC_THEORY_NEXT_CAMPAIGN_20261014.md).
TH-001–005 are proposed work, not an executed campaign. Empirical/source
gates and the original cost/model distinctions remain unchanged.

Bounded diagnostics:

```bash
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 \
python3 fin_theory_synthesis_20261014/verify_bridges.py \
  --output /tmp/fin-theory-check.json
```

Run from repository root. It reads sources and writes only the chosen output
and its accompanying index. Optional --tool-replay points to a READY scratch
replay. Use the source's unit tests only on a copy because they write fixtures
and output records. Numerical checks do not prove the future stochastic limit.
