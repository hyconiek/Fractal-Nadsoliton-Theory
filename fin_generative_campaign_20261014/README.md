# FIN RG-C01–RG-C05 — executable research handoff

**14 October 2026.** Locally executed bounded finite-N q12 FIN study, based on the supplied plan and public versioned formulas. No repo checkout was available; see `SOURCE_ACCESS.md`.

## What is inside

- `HANDOFF.md`, `CLAIM_REGISTER.json`, `EVIDENCE_LEDGER.md`, `SOURCE_ACCESS.md`.
- `SOURCE_PLAN_20261014.md`: the original user-provided plan, unchanged.
- `RG-C01` to `RG-C05`: reports, structured results, nonconclusions, replay and hashes, plus task-specific artifacts.
- `RG-C02/INTERVAL_CERTIFICATE.json` and `INTERVAL_PROOF.md`: exact-decimal-kernel, finite-alphabet outward enclosure **conditional on mpmath.iv**.
- `RG-C03/PROCESS_MAP.json`, `MEMORY_REPORT.json`: source generator and finite-memory comparisons.
- `RG-C04/COUNTERMODEL_COMPARISON.csv`: normalized countermodel controls.
- `RG-C05/PREREGISTRATION.json`, `PREDICTIONS.json`, `VALIDATION.json`: held-out synthetic model consistency.
- `run_campaign.py`, `interval_certify.py`, `independent_verify.py`: reproducible producer, independent arithmetic certificate, independent second verifier.
- `INDEPENDENT_VERIFICATION.json`: 14 separately implemented checks, all passed.
- `RG-C06_CHECKPOINT_FOR_ARCHITECT.md`: decision **not made** by executor; `RG-M00_OPEN_INPUT.md` identifies unavailable inputs.

## Minimal reproduction

Install `numpy`, `scipy`, `mpmath` in a Python 3 environment and run from this directory:

```bash
OPENBLAS_NUM_THREADS=1 python run_campaign.py
OPENBLAS_NUM_THREADS=1 python interval_certify.py
OPENBLAS_NUM_THREADS=1 python independent_verify.py
sha256sum -c MANIFEST.sha256
```

**Note:** the first three commands regenerate output files and may alter the checksums in this original manifest because runtime values are recorded; run `sha256sum -c MANIFEST.sha256` **before** replay to verify the distributed archive. After replay regenerate the manifest with the included `finalize_package.py` only after source plan availability is confirmed. `finalize_package.py` requires the original uploaded plan path and is provided for traceability, not needed for running numerical checks.

For manual verification, start with `HANDOFF.md`, `RG-C02/INTERVAL_PROOF.md` and `INDEPENDENT_VERIFICATION.json`.
