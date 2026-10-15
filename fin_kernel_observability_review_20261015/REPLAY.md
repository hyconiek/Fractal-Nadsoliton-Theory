# KO intake replay

From the repository root with numpy, scipy and mpmath available:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 fin_kernel_observability_review_20261015/certify_event.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 fin_kernel_observability_review_20261015/verify_intake.py --replay --output /tmp/fin-ko-reverification.json
```

The certificate uses accepted rational spectral boxes, an independent input
DFT, positive-series erf bounds and exact raw-dyadic outward rounding.
It assumes mpmath.iv correctness, not a proof-assistant backend.

The verifier checks coverage including the producer's separately excluded
PACKAGE_VERIFICATION.json, two matching archive wrappers, real prior reviews
and saved inputs. It copies sources to a cache-free temporary directory,
runs the source package checker before scientific writers, and rechecks
all original bytes afterwards. Source numerical output is not cross-environment
byte-identical. build_package.py/finalize_package.py are not replay requirements
and must not rewrite the historical source.

After intentional edits, document/package validation is:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 fin_kernel_observability_review_20261015/finalize_review.py
```

This updates only this review's packaging, not scientific source or DKO results.
Read PROOF_COMPLETIONS_AND_ARITHMETIC.md separately: finite assertions and
hashes do not mechanically certify the asymptotic CLT/local-window argument.
