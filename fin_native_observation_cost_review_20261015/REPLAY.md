# NOC intake replay

From the repository root with numpy/scipy/mpmath:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 fin_native_observation_cost_review_20261015/certify_mean_gap.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 fin_native_observation_cost_review_20261015/verify_intake.py --replay --output /tmp/fin-noc-reverification.json
```

The verifier binds incoming source coverage, mapped actual external inputs,
prior review manifests and scientific inputs. It runs the producer package
checker before its scientific writer on a cache-free temporary copy, then
checks all original bytes unchanged. finalize_noc.py is not a replay
requirement: it authors reports and uses historical absolute paths.

Additional checks use actual finite labelled jumps, exact private sampling
moments, full q12 N2 Duhamel quadrature and integer-root floor identities.
These are bounded diagnostics, not asymptotic proofs or tiny-domain trapping
experiments. The gap/tail certificate has explicit spectral and iv assumptions.

After intentional review edits, validate/package only this review:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 fin_native_observation_cost_review_20261015/finalize_review.py
```

No SRN task, large-N trajectory or physical experiment is run by these checks.
