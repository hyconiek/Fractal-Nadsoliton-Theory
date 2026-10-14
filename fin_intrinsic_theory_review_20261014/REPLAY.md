# Replay and proof inspection

From the repository root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 fin_intrinsic_theory_review_20261014/verify_intake.py --replay --output /tmp/fin-th-intake-verification.json
```

The verifier copies incoming writers to a fresh `/tmp` directory, preserves
archived source files and invokes four small scripts. `write_reports.py` is
not run. It independently binds versioned model formulas, the local fold and
spectral enclosures, checks the full tangent spectrum and replays R7P-031's
small interval calculation. The optional historical `.pyc` entry is not needed
to execute science and may be absent in a versioned checkout.

Required Python packages: numpy, scipy, mpmath. Different versions can change
byte hashes; record numerical drift instead of promising byte identity.

The fresh literature capture `/tmp/fin-th-intake-kraaij.md` is optional for
replay. Its historical hash is evidence of this intake's access, not a required
reconstruction of the author's web browsing. Inspect the full theorem at
[the primary source](https://arxiv.org/html/1502.05185v3#S2.Thmtheorem3)
against `INTAKE_20261014.md` and the incoming TH-002 proof.

Manually inspect the analytical flux-duality, quasipotential connection and
Fourier-recovery arguments; no finite check proves them by itself. Inspect
`FAST_MODE_PROOF_ROUTE.md` only as a proof plan. Source and intake PASS flags
do not prove the stopped stochastic-fold convergence target.
