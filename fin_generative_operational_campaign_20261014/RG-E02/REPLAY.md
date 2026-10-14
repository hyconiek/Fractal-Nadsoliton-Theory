# Replay

From the campaign root with Python 3, NumPy, SciPy, mpmath available:

```bash
OPENBLAS_NUM_THREADS=1 python run_operational.py
OPENBLAS_NUM_THREADS=1 python interval_static.py
OPENBLAS_NUM_THREADS=1 python independent_verify.py
```

`run_operational.py` regenerates prospective model probabilities and diagnostics; `interval_static.py` separately reconstructs the supplied strict-kernel static probabilities in outward intervals (conditional on mpmath.iv); `independent_verify.py` reconstructs FIN with an independent X X^T basis and dense matrix exponentials. The `make_reports.py` script reconstructs report and manifest metadata, so do **not** run it before verifying the distributed source MANIFEST. All reference data are mathematical/synthetic.
