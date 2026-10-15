# Clean scientific replay (new RI work only)

From this directory with Python 3, NumPy and SciPy installed:

```bash
python independent_checks.py
OPENBLAS_NUM_THREADS=1 python deterministic_table.py 0.10 0.10 0.025 4.5 3.8 > /tmp/GRID_DIAGNOSTIC_COARSE.json
OPENBLAS_NUM_THREADS=1 python deterministic_table.py 0.075 0.075 0.0125 4.5 3.75 > /tmp/GRID_DIAGNOSTIC_FINE.json
OPENBLAS_NUM_THREADS=1 python deterministic_table.py 0.06 0.06 0.0125 4.5 3.78 > /tmp/GRID_DIAGNOSTIC_REFINED.json
python verify_package.py
```

Compare reconstructed JSON numerically with `GRID_DIAGNOSTIC_*.json` for sensitivity; cross-version floating byte identity is **not** a theorem. All numerical grids are uncertified approximations. `independent_checks.py` writes `INDEPENDENT_CHECKS.json` and exact Fraction tests. `verify_package.py` validates every file against `MANIFEST.sha256`, excluding the manifest itself. Previous campaign sources and code are NOT embedded; references and actual accessible hashes appear in `SOURCE_HASHES.json`.
