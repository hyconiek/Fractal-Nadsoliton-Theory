# Reproduction sequence

1. Inspect `SOURCE_ACCESS.md`; understand that commit provenance is incomplete.
2. Verify package `MANIFEST.sha256` with `sha256sum -c MANIFEST.sha256` in the shipped state (before executing code).
3. Install Python 3 + numpy/scipy/mpmath, run `OPENBLAS_NUM_THREADS=1 python run_campaign.py`.
4. Run `OPENBLAS_NUM_THREADS=1 python interval_certify.py`; do not use an unqualified binary64 bound for exact-kernel formal scope.
5. Run `OPENBLAS_NUM_THREADS=1 python independent_verify.py`, require `passed=true`.
6. Review every `RG-C0x/RESULTS.json`, `NONCONCLUSIONS.md`, and `INPUTS_UNAVAILABLE.md`.
7. If a full repository clone becomes available, bind actual latest Git SHA and original source SHA, inspect diffs, and rerun before attaching findings to upstream.
