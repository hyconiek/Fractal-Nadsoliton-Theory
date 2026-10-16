# Replay of bounded controls

From repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 fin_balanced_phase_cut_research_20261016/finalize_research.py --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 fin_balanced_phase_cut_research_20261016/verify_research.py --output /tmp/fin-balanced-phase-replay.json
```

The verifier uses stdlib exact fractions and integer determinants for its
finite controls, and separately numpy/scipy through versioned safe_seed for
nominal A7 rate normalization. A separate q2 binomial check tests soft-mesh
conventions, not q12 critical convergence. No large lattice is enumerated.
Proofs of actual growing-N statements are in the four theorem files, not
established by these bounded controls. Changing floating last digits across
environments is not a mathematical failure of an analytical theorem.

`finalize_research.py` without --check writes metadata/hashes only for this
new package after validating the existing results. It does not regenerate
or alter any predecessor/source package or count as scientific evidence.
