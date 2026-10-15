# DKO intake replay

From repository root with numpy/scipy/mpmath:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 fin_dynamic_kernel_observability_review_20261015/certify_sign_event.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 fin_dynamic_kernel_observability_review_20261015/verify_intake.py --replay --output /tmp/fin-dko-reverification.json
```

The manifest contains a historical .pyc. The package checker therefore runs
on an **archival integrity-only copy**, while scientific independent_verify.py
runs on a distinct cache-free copy. The historical bytecode is never imported
or executed. If a clean checkout omits that cache, the intake reports the
producer's cache-only package failure separately; it does not reconstruct
scientific evidence from bytecode or waive any missing scientific file/hash.

All original sources, including separately excluded PACKAGE_VERIFICATION.json,
are checked unchanged afterwards. Numerical JSON need not be byte-identical
across environments. Replayed flags are bounded diagnostics, not a mechanical
proof of the stochastic or Gaussian information theorems.

The new event certificate uses rational spectral inputs and a positive arcsin
series with geometric tail and exact outward decimal rounding. Its backend
assumptions and nonconclusions are explicit. It does not give finite-N errors.

After intentional review edits, validate/package this review only:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 fin_dynamic_kernel_observability_review_20261015/finalize_review.py
```

No NOC task, large-N trajectory or external experiment is run by these checks.
