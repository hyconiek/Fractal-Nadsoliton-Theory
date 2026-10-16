# Reproducibility

Run `PYTHONDONTWRITEBYTECODE=1 python3 verify_kr.py` from this directory. It regenerates `FINITE_EXACT_CHECKS.json` deterministically using Python's standard-library `fractions` and validates finite q4/N2 rate algebra. Run `PYTHONDONTWRITEBYTECODE=1 python3 verify_package.py` from the same directory to check the SHA-256 manifest and generated-file coverage. This replay does *not* claim to reproduce the UF/MT asymptotic proofs, import a current git checkout, or reconstruct unavailable predecessor bytes.

The q4 wheel is a deliberately small **separate algebra test**, not a different asymptotic q12 critical diffusion.
