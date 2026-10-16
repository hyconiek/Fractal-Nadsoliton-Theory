# Reproduce the HM intake

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 fin_auxiliary_mediator_review_20261016/verify_intake.py --replay --output /tmp/fin-hm-intake-verification.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_auxiliary_mediator_review_20261016/finalize_review.py --check
```

Requires Python3, numpy, scipy and sympy. Source manifest coverage is checked
before copying; package_check.py runs before verify_hm.py on a temporary
cache-free copy. The source verification writer never touches archived output.
Do not run finalize_hm.py as replay: it authors status/narrative/input files
and depends on original unavailable upload paths.

The separate verifier binds actual MP7/model/MT/NT/KI/DKO files and uses fresh
finite algebra, the full18D local Jacobian, nonzero exact compression and
234-state reset-lift moment/intertwining fixtures. Uniformization operations
are floating, not an interval probability certificate. The three residual
atoms test reset algebra, not the full Gaussian conditional law.

Stationary critical convergence is an analytical theorem in the proof
completion, not a consequence of finite tests. Source descriptive flags are
not independently paid by rerunning the same True expression.

finalize_review.py without --check refreshes only this review's mechanical
metadata/manifest; --check is read-only. No old source/report is rewritten.
