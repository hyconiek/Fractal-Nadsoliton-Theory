# Bounded KI intake replay

From the repo root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 fin_kinetic_invariants_review_20261014/verify_intake.py --replay --output /tmp/fin-ki-intake.json
```

Requirements: numpy and scipy. The verifier checks the original72-entry
manifest,8 source inputs,37 task inputs and61 SF archive members; it reuses
unchanged local SF inputs/completion. Writers run in a fresh `/tmp` copy.
Package hash checks precede source writers, not follow modified numerical files.
Authored build_reports/finalize scripts are never run.
The replay copy omits historical `__pycache__` directories; archived source
files are not deleted or changed.

The independent intake checks actual rate derivatives, strict tangent inertia,
mobility/noise, nonlinear Onsager identity,396 finite DB cases and the actual
quadratic-square generator identity. It uses versioned safe_seed and different
QR/Kronecker algebra, not the incoming producer's numerical core.

Read KI-001/002/003 proofs and the locally available SF joint-clock completion
separately. Numerical replay does not prove stochastic moment/weak-limit bounds.
Source flag defects and the thermal-native-stop clarification are in INTAKE.

Uniform critical spectral bounds are exact Fraction propagation of previously
admitted boxes. No uniform critical path campaign is executed by this verifier.
