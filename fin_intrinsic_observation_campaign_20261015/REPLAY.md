# Computational replay and dependency provenance

From a Python environment with numpy/scipy:

```bash
python run_observation.py
python independent_verify.py
python verify_package.py
```

All commands run in this package directory. The `source_context/theory_core.py` is a numerical source reimplementation copied byte-for-byte from the earlier MT ZIP, not an interval oracle or independent A7 derivation. The verifier does **not** import main observation code nor that theory-core; it independently reconstructs A from the exact-decimal parameter strings in double precision, and uses Python `Fraction` for small exact hypergeometric probability checks. All double-precision assertions have explicit tolerances.

The archived `source_context/FIN_MT_FULL_HANDOFF.zip` is unchanged, and the user-pasted OBS plan/MT intake are stored with original hashes, not normalized into historical source files. The accepted MT reviewer `PROOF_COMPLETIONS.md` was not provided as full text, so the current campaign treats its described completion as an inherited premise rather than claiming to re-audit it.

No source AGENTS, Git checkout or remote repository state was modified. Numerical replay hashes may vary across scipy/numpy versions; equations, exact rational checks and disclosed error tolerances are the primary reproducibility targets.
