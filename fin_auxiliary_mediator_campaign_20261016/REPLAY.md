# Reproduction

From this directory: `OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 verify_hm.py`. For scientific verification, audit the exact proof sections in HM_THEOREM_PROOFS.md and compare the boundary/moment ledger with the actual current MT/NT reviews. Manifest: `sha256sum -c MANIFEST.sha256` from this directory. `package_check.py` also validates the expected manifest members. No past package is required to run small numerical checks, but independent source custody is explicitly limited.
