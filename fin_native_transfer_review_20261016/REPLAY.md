# Replay the NT intake

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 fin_native_transfer_review_20261016/verify_intake.py --replay --output /tmp/fin-nt-intake-verification.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_native_transfer_review_20261016/finalize_review.py --check
```

Python3, numpy, scipy and sympy are required. The intake verifies incoming
manifest coverage and separately binds the excluded package checker/metadata.
Source package_check.py runs on a copy before verify_nt.py, which writes
only that copy's verification.json. No source archive or old review is changed.

The separate verifier reconstructs the labelled parent from versioned safe_seed,
uses orbit intertwining to obtain the count generator and propagates by a
Poisson uniformization series. Floating error/tail is not an all-operation
interval probability certificate. Complex chronological cylinders and symbolic
checks validate conventions, not a stochastic asymptotic theorem by sampling.

Actual hashes bind the MT/OBS stationary proof, prior research and model.
HM cavity/mobility checks are bounded seeds for the new plan, not execution
of HM or physical adoption of its diagnostic law.

finalize_review.py without --check refreshes only this review's mechanical
metadata/manifest after an authorized edit; --check is read-only. The old
source's missing uploaded inputs are not reconstructed or marked matched.
