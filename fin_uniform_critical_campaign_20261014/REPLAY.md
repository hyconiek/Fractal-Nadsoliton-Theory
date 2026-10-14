# Reproduction steps

1. Check `MANIFEST.sha256` before running any writer; it hashes the original frozen reports and all bundled inputs.
2. Create an environment with Python 3, numpy and scipy (see ENVIRONMENT.json).
3. Execute `python uf_analysis.py`, then `python independent_verify.py`. Compare numerical values to DIAGNOSTICS.json and INDEPENDENT_VERIFICATION.json within environment-appropriate floating tolerances; exact binary64 hash replay is not promised.
4. The independent verifier reconstructs A from circulant spectral projectors without importing the producer's source model; both use the same mathematical definitions.
5. Review UF-001 to UF-004 proof documents for analytical obligations. Script checks do not prove stopped process convergence.
6. Do not automatically rerun `make_reports.py`, which re-authors verdicts and hashes. The previous KI ZIP and original user plan/intake are preserved unchanged.

No live remote repository access, external experiment, or Git mutation occurred.
