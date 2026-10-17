# NS intake replay

From the repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 fin_nonlinear_stability_transfer_review_20261017/verify_intake.py --replay --output /tmp/fin-ns-intake-replay.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_nonlinear_stability_transfer_review_20261017/finalize_review.py --check
PYTHONDONTWRITEBYTECODE=1 python3 fin_coupled_evolution_source_review_20261017/finalize_review.py --check
git diff --check
```

The producer REPLAY/verify.py writes exact_checks.json and must run only on
a fresh temporary copy. The intake does this, checking source coverage before
and after. Recorded byte equality applies to this replay, not all platforms.
The independent controls use different source-type matrices and rational
signed-wall bounds, with the accepted safe polynomial parser as a tool.

Analytical min-max, energy/continuation and transfer estimates are written
separately; no finite fixture proves them. The finalizer writes only this
new review's own index/manifest in write mode; --check is read-only.
Current source hashes must match. Missing original wrapper is not a mismatch
and is not recreated as historical custody. No browser/physical experiment.
