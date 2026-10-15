# Replay of the OBS intake

Run from the repo root, with numpy, scipy and mpmath available:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 fin_intrinsic_observation_review_20261015/verify_intake.py --replay --output /tmp/fin-obs-reverification.json
```

The checker verifies incoming hashes and actual predecessor dependencies,
creates a fresh cache-free `/tmp/fin-obs-intake-*/campaign` copy, runs the
source package checker before scientific writers, and rechecks original
sources afterwards. The expected prototype-support exception is recorded
as an isolated defect; the corrected helper is tested independently.

Do not run `build_metadata.py` or `finalize_package.py` against the original
campaign. They author reports/manifests and do not add scientific evidence.

Package/document validation, after any intentional review edit:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 fin_intrinsic_observation_review_20261015/finalize_review.py
```

This regenerates only this review's packaging objects. It is not an OBS or
RI theorem execution. To retain the original intake evidence, direct fresh
scientific output to `/tmp` and compare rather than overwrite verification.json.

Review the proofs and [corrections](PROOF_COMPLETIONS_AND_CORRECTIONS.md)
separately: bounded fixtures are not a stochastic convergence certificate.
The Whitt capture is an ignored local web cache; the public reference and
its accessed-text hash are retained even if that cache is absent elsewhere.
