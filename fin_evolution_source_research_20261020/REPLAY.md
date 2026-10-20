# Bounded replay

From repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 timeout 60s python3 fin_evolution_source_research_20261020/verify_research.py --output /tmp/fin-source-probe-replay-20261020.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_evolution_source_research_20261020/check_package.py --check
PYTHONDONTWRITEBYTECODE=1 python3 fin_complete_anonymous_history_review_20261020/finalize_review.py --check
PYTHONDONTWRITEBYTECODE=1 python3 fin_post_rp004_history_review_20261020/finalize_review.py --check
PYTHONDONTWRITEBYTECODE=1 python3 fin_strategy_audit_20261019/verify_strategy.py --inputs-only
```

The new verifier uses Fraction arithmetic and independently constructed
finite block laws/windows. It does not run old producer writers or simulate
long histories. The native p check imports versioned safe_seed and is
explicitly binary64/nominal, not a transcendental input enclosure.

The full conditional proofs are in PROOF.md. The number of fixture checks
does not prove stationarity, the all-order statement or physical source
closure. External literature is context only, not a replay dependency.

`check_package.py` writes this new package's input binding and manifest only
without `--check`. With `--check` it is read-only. It preserves every old
source/review, historical strategy and AGENTS. Its success is not success
at finding F01: the claim register must retain OPEN_SOURCE_NOT_DERIVED.
