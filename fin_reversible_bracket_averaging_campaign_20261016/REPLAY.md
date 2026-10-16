# BA reproduction and boundaries

From this directory (with Python and SymPy available):

```
PYTHONDONTWRITEBYTECODE=1 python3 verify_ba.py
PYTHONDONTWRITEBYTECODE=1 python3 verify_package.py
```

`verify_ba.py` independently builds a rational reversible six-state ring (not q12), reproduces flux constraints and mean brackets, solves exact Poisson equations for several positive `tau`, and validates the exact reduced 3-component leading coefficient. It records `EXACT_DIAGNOSTICS.json`; if that file is regenerated, its byte identity is environment-dependent only in formatting and deterministic rational algebra, not a proof of an asymptotic q12 family.

`verify_package.py` checks retained `MANIFEST.sha256` and basic JSON source/claim/fixture integrity. The manifest excludes itself and `PACKAGE_VERIFICATION.json`; the latter is generated after the manifest and hash-bound separately by the archive wrapper check. The package contains no cached native predecessor files or automatic AGENTS edits.

The main results (cut flux and component formula; Poisson maximal estimate) are analytic and should be independently audited using the actual current KR review. The old AGENTS.md available locally is a historical snapshot, not verified live HEAD. Missing original sources have no claimed custody.
