# Minimal replay and verification

1. Verify distributed artifact integrity **before** running any writer: `sha256sum -c MANIFEST.sha256` from this directory (Linux/macOS).
2. Python 3 with NumPy/SciPy/mpmath:

```bash
OPENBLAS_NUM_THREADS=1 python run_operational.py
OPENBLAS_NUM_THREADS=1 python interval_static.py
OPENBLAS_NUM_THREADS=1 python independent_verify.py
```

3. Inspect JSON under RG-E02, RG-E03, RG-E04; verify 58/58 checks and conditional N2 static enclosures.
4. Regenerate static narrative/manifest only if desired: `python make_reports.py`. Replays may change floating serialization across OS/numpy/SciPy; compare model laws with numeric tolerances, not byte identity.
5. For independent actual experiment, do not execute without separate authorization and calibration. This pack has no physical readout implementation.
