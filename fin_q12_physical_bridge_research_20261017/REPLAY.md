# Replaying physical-compatibility controls

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 fin_q12_physical_bridge_research_20261017/verify_bridges.py --output /tmp/fin-q12-bridge-replay.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_q12_physical_bridge_research_20261017/finalize_research.py --check
```

The checks use SymPy exact algebra, rational finite-bank/partial-SWAP
controls and versioned safe_seed numerical diagnostics via numpy/scipy.
They do not prove a physical model or mechanize the analytical statements.
The browser capture log binds three actual successful public-page reads.
Optional ignored cache is checked if present; its absence on another checkout
is reported rather than reconstructed. Replaying the verifier does not
pretend to repeat the historical browser visit.

browser_capture.py is an explicit helper for browser-harness globals.
Only task-created public tabs on an allowlist are captured; no cookies,
credentials, browser profiles, recordings or unrelated tabs are exported.
The package finalizer writes only this new package's metadata/manifests.
