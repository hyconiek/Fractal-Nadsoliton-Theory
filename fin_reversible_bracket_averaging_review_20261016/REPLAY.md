# Source-preserving BA intake replay

From repository root, with SymPy available for the incoming producer only:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 fin_reversible_bracket_averaging_review_20261016/verify_intake.py --replay --output /tmp/fin-ba-intake-replay.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_reversible_bracket_averaging_review_20261016/finalize_review.py --check
```

The intake checker uses stdlib arithmetic and imports only the older
versioned AC finite-model helper, never incoming BA code. It checks current
hashes, archived members, exact nonuniform count cuts/Dirichlet pairings,
cycle-current Poisson costs/residue and seed symmetry/dual fixtures.

With --replay it copies the source to /tmp, runs verify_package.py before
verify_ba.py, and compares the regenerated exact-data bytes. Original files
are checked again after replay. PACKAGE_VERIFICATION.json is producer-excluded
and bound separately, not silently counted inside the52 entries.

Historical AGENTS bytes remain unavailable. Actual current dependencies and
older saved append-only ledger inputs are checked separately by hash/prefix.
No missing custody or uploaded-plan identity is reconstructed.

The finalizer is mechanical. --check is read-only. Without it only this
review's PACKAGE_CHECKS.json and MANIFEST.sha256 are generated. It does not
author scientific narratives or modify prior sources.
