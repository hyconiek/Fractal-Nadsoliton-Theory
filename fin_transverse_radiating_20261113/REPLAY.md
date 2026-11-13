# Replay

From the repository root with Python and NumPy:

    OPENBLAS_NUM_THREADS=1 python fin_transverse_radiating_20261113/run.py

This overwrites only this package's four CSV traces and results.json.
The parent static package supplies Delta_6. No historical generator is run.
NumPy FFT and RK4 are shared methods with prior work; this is not an
independent-team replication. INPUTS.json and MANIFEST.json record hashes.
