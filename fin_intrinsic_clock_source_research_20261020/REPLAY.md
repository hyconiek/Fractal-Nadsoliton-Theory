# Replay and integrity

From repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 fin_intrinsic_clock_source_research_20261020/verify_research.py --check --output /tmp/fin-intrinsic-clock-check.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_strategy_audit_20261019/verify_strategy.py --inputs-only
git diff --check
```

The new verifier uses stdlib exact Fraction arithmetic, formal prime-factor
logarithms and an independently implemented finite ES lift. It solves actual
killed fibres rather than merely comparing the displayed exit formula with
itself. q12 checks enumerate 1728 block words, not a count lattice or large
hidden-state graph. No numerical quadrature proves the infinite series.

The 332 controls include repeated per-pair/per-word assertions; they are not
332 independent proofs. Source/AGENTS input hashes and predecessor manifests
are checked separately. If later instructions are appended, preserve this
historical full-file input and use an explicitly labelled historical context
for its replay; never rewrite saved scientific hashes or follow the archived
AGENTS as live authority.

Authoring --record writes only this fresh directory's derived JSON, input
ledger and manifest. Subsequent --check writes only the requested /tmp
result. No old checker, writer, narrative or manifest is altered.
