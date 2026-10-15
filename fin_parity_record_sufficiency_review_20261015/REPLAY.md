# LI intake replay

From the repository root with numpy/scipy:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 fin_parity_record_sufficiency_review_20261015/verify_intake.py --replay --output /tmp/fin-li-reverification.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_parity_record_sufficiency_review_20261015/kernel_family_certificate.py
```

The verifier checks source coverage, resolves named external references and
reports the unavailable uploaded plan separately. It checks actual prior
review manifests and saved inputs, then performs bounded q12/exact toy checks.
The producer checker runs before its writer on a cache-free temporary copy.
The source package is verified again afterwards. The root producer REPLAY.md
is absent; its per-task notes remain preserved.

The source's20527 assertions are floating q3 diagnostics, not exact rational
or A7 process-limit certificates. The proof audit in PROOF_COMPLETIONS.md is
a separate analytical assessment. No large simulation or prior campaign is rerun.

After intentional review edits, validate/package this review only:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 fin_parity_record_sufficiency_review_20261015/finalize_review.py
```

This is document/hash validation, not mechanized LI proof or KO execution.
