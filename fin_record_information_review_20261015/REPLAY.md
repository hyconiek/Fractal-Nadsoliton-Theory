# RI intake replay

From the repo root with numpy and scipy:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 fin_record_information_review_20261015/verify_intake.py --replay --output /tmp/fin-ri-reverification.json
```

The checker verifies incoming hashes and local predecessor dependencies,
reconstructs small actual q12 generators and likelihoods independently, and
creates a cache-free temporary source copy. There it runs the source manifest
checker **before** writers, the 5626 rational assertions, and all three
deterministic grids. Source files are checked again afterwards.

Do not run `assemble_package.py` against the original source or treat replay
grid differences as rigorous error bounds. All finite numerical fixtures are
algebra diagnostics, not large-N trajectories in the tiny MT domain.

Document/package validation after intentional review edits:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 fin_record_information_review_20261015/finalize_review.py
```

This writes only the review packaging objects. It neither proves the
analytical arguments nor executes LI. The primary filtering reference URL
and accessed-text hash are retained in verification.json; the web capture
is an ignored local cache, not mandatory historical producer evidence.
