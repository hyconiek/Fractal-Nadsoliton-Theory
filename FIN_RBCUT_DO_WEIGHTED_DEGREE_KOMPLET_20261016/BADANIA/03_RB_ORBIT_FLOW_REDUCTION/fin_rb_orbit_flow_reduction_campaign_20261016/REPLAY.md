# Replay and custody

Run in a clean Python >=3.10 environment (stdlib only):

```bash
PYTHONDONTWRITEBYTECODE=1 python3 verify_orbit.py
```

The script uses exact `Fraction` arithmetic for the 36-node and 6-node toy fixtures, generates `EXACT_CHECKS.json` deterministically and fails on any violated stated equality. It explicitly checks all 12 dihedral actions on the toy, all 64 quotient cuts, a primal lifting, a dual rate-cap certificate, a modified-graph infeasibility control, and independent exhaustive 64 full / 8 orbit cuts in a smaller example.

To verify the archive, extract it safely to a temporary directory and compare every SHA-256 listed in `MANIFEST.sha256` from the package root. `INPUTS.sha256` lists hashes of source files as available in the current container, not the user's latest Git HEAD. The archive is not a source checkout or a peer-reviewed formalisation.
