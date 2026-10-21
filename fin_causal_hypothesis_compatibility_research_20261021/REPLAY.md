# Read-only replay

From the repository root:

```sh
prlimit --as=1073741824 --cpu=55 -- timeout 55s env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 fin_causal_hypothesis_compatibility_research_20261021/verify_research.py --check
git diff --check
```

The checker recomputes finite diagnostics, validates its exact coefficients,
compares actual input hashes and predecessor manifest members before/after,
and checks all own non-cache members and coverage. No old writer/finalizer
is executed. It does not assert that old full-file AGENTS input hashes match
today's live ledger. The separate strategy inputs-only command was actually
run and reported 87 bound files with its new-live-context notice.

`--diagnose` performs the same bounded scientific checks without authoring
metadata. `--record` is only for authoring this NEW packet's own RESULTS,
INPUTS, scoped duplicate-search output and manifest, after source-preservation
checks. Do not use it to repair a changed source input or to reseal history.
Run record and check sequentially to completion, not as concurrent writers.

Thirty-one exact Fraction/structural controls plus three ordinary binary64
diagnostics pass. The largest exact linear solve is 27 states. One nonlinear
q3/N3 fixture at +/-1e-4 checks a centred derivative, not a parameter sweep
or an outward enclosure. Versioned safe_seed supplies only a nominal A7
arithmetic check; the actual positive sign in the theorem uses PSD and the
parity eigenvalue analytically. No browser/literature cache is required.
