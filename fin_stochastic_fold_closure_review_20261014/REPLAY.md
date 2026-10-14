# Bounded verification and proof reading

From the repo root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 fin_stochastic_fold_closure_review_20261014/verify_intake.py --replay --output /tmp/fin-sf-intake.json
```

Packages: numpy/scipy for diagnostics. mpmath is additionally needed for the
source's quadrature cross-check. Exact public exit-endpoint verification is
stdlib Fraction arithmetic. Scripts write their outputs only in a temporary
campaign copy. Source manifest checks run before writers; no authored metadata
builder is run, and the original60-entry manifest is rechecked afterwards.

The verifier reuses unchanged predecessor input hashes and actual fold bytes,
rather than redoing a root census. It independently checks the chart/basis,
Lyapunov matrix through a Kronecker solve and actual finite-jump quadratic terms.
It rejects a deliberately inward exact endpoint. Source41/41 and intake26/26
are computational flags, not mechanically checked analytical theorems.

Manually read SF-001 FAST_CHART_PROOF, SF-002 GENERATOR_ENERGY_PROOF,
SF-003 STOPPED_CONVERGENCE_THEOREM, and the intake's
STOPPED_LIMIT_COMPLETION. Admission includes that completion; source interior
tests alone do not license absorption. The primary UI result is
[Whitt Theorem5.3](https://arxiv.org/html/0712.1929v2#S5.Thmtheorem3),
not a blind deterministic-clock FCLT at a random stop.

The optional `/tmp/fin-sf-martingale-primary.md` capture is access metadata,
not a required artifact or proof of the producer's historical browsing.

KI execution is separate and has not started.
