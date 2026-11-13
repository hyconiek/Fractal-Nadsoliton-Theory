# Replay

From the repository root with Python, NumPy and SymPy:

    python fin_transverse_kink_stability_20261113/check.py

This recomputes results.json. The checker verifies the exact spatial
factorization and eigenfunction, compares the 12 by 12 matrix spectrum with
the Fourier calculation, checks multiplicities and the numerical gap.
It does not execute historical source generators or establish interval bounds.
INPUTS.json pins historical inputs; MANIFEST.json pins this package's files.
