# Independent replay and integrity

Run from extracted directory (no network necessary):

```bash
PYTHONDONTWRITEBYTECODE=1 python3 verify_rb.py
PYTHONDONTWRITEBYTECODE=1 python3 verify_package.py
```

`verify_rb.py` checks integer Fourier/parity conventions for the q12 label map and an exact synthetic 12-vertex rational native-flow/dual fixture. Its result is bounded finite diagnostics and does not establish any true q12 large-N feasibility theorem.

`verify_package.py` verifies manifest SHA-256 coverage of all RB artifacts other than self-referential metadata. Source scientific documents are *not* silently included; to audit them, retrieve the exact upstream bytes and compare with the digests in `SOURCE_PROVENANCE.json`/`INPUTS.sha256`. Missing newer review sources cannot be proved identical by this package.

The only accepted model-specific new formula is an exact finite identity conditional on the already specified native generator and support. The general directed graph maxflow theorem is proved in RB-002. The unresolved all-N native ideal-sign property remains `OPEN`.
