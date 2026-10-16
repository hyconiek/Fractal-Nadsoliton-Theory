# Replay

Run from a clean copy of this new directory:

`python3 verify_weighted_degree.py`

then

`python3 verify_package.py`

Both require only Python standard library. The first recomputes exact Fraction toy and small-network boundary fixtures and rewrites `FINITE_CHECKS.json`. The second recomputes hashes of all files listed in `MANIFEST.sha256`. Hashing and exact toy flags test reproducibility, **not** a formal proof of the analytic all-N FIN theorem. The full proof is readable in `UNIFORM_NATIVE_EDGE_ELLIPTICITY.md` and `WEIGHTED_ORBIT_CUT_EQUIVALENCE_THEOREM.md`.

`INPUTS.sha256` binds actually accessed source files and archived MT/UF wrappers; inaccessible current review objects have no fabricated hash. No upstream authorial builder is run, no historical source changed, and no external data are used.
