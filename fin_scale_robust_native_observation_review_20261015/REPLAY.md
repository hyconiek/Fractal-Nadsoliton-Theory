# Reproduce the SRN intake

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 fin_scale_robust_native_observation_review_20261015/certify_angular_gap.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 fin_scale_robust_native_observation_review_20261015/verify_intake.py --replay --output /tmp/fin-srn-verification.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_scale_robust_native_observation_review_20261015/finalize_review.py --check
```

Requires Python3, numpy, scipy and mpmath; the copied producer also uses sympy.
The intake verifies the incoming manifest before copying; only verify_srn.py
runs on a cache-free temporary copy and writes its local verification output.
Do not run finalize_srn.py against the historical source; it rewrites ledgers
and expects original unavailable uploads.

The saved verification binds original source bytes, actual predecessor inputs,
logs, replay byte status, exact finite fixtures and the ideal arithmetic
certificate. The new scientific proof is analytical in PROOF_COMPLETIONS.md;
small checks do not prove its all-N hypotheses. No massive simulation is needed.

finalize_review.py without --check refreshes only this review's mechanical
package metadata/manifest after an authorized review edit. --check is read-only.
