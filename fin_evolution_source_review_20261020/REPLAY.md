# Replay and source preservation

From repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 fin_evolution_source_review_20261020/verify_intake.py --replay --output /tmp/fin-causal-source-final-20261020.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_evolution_source_review_20261020/finalize_review.py --check
PYTHONDONTWRITEBYTECODE=1 python3 fin_evolution_source_research_20261020/check_package.py --check
PYTHONDONTWRITEBYTECODE=1 python3 fin_complete_anonymous_history_review_20261020/finalize_review.py --check
PYTHONDONTWRITEBYTECODE=1 python3 fin_strategy_audit_20261019/verify_strategy.py --inputs-only
```

The new verifier checks complete ES source coverage and its fourteen inputs.
It copies ES into a temporary sibling-input tree, symlinking existing
dependencies; the source checker is invoked read-only **before** the writer.
Only the copied RESULTS is regenerated, with 60-second/1-GiB child limits.
No old producer writer is run against the repository. The full old source
snapshot and live AGENTS are rechecked after replay.

New exact rational checks use an independently built sparse 81-state
analysis source, not ES scientific imports. Filter polynomials, impossible
events, departure survival, source/visible defect distinctions and entropy
inequalities are bounded fixtures, not mechanized proofs. The separately
rebuilt 144-state native parent uses ordinary floating arithmetic and is
not an outward certificate. No long stochastic simulation or network access
is required for this continuation.

Without flags, the new finalizer seals **only this new review** and its
root status documents. With `--check` it is read-only. Earlier source,
review, strategy, AGENTS and their historical hashes are never rewritten.
