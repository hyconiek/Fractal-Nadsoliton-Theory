# Reproducibility

From the new handoff folder:

`PYTHONDONTWRITEBYTECODE=1 python3 verify_cut_shell.py`
`PYTHONDONTWRITEBYTECODE=1 python3 verify_package.py`

`verify_cut_shell.py` rebuilds independent exact q12 Fourier/one-step identities and an exact small transshipment negative control; no random seeds, nonstandard packages, or producer imports. Its PASS flags are bounded finite **diagnostics**, not all-N mathematical proof.

The mathematical proofs are written explicitly in `NATIVE_SHELL_THEOREM.md`. The moving-shell stationary statement additionally assumes the previously admitted MT-003 local Gibbs weak convergence; its original archived member hash is included in `INPUTS.sha256`. Do not substitute original review statements that were not locally available.
