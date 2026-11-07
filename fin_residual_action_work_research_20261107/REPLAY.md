# Replay

From the repository root, with Python 3, NumPy and SymPy installed:

    python fin_residual_action_work_research_20261107/check_work.py

The script checks the pinned predecessor hashes and writes verification.json.
The recorded environment used NumPy 2.3.5 and SymPy 1.14.0. Floating-point
results may differ in the last digits across FFT implementations.

The base grid has 4096 points on [-256,256). Refinements use 8192 points
on that interval and 6144 points on [-384,384). Times are integrated
analytically, not by time stepping. The spatial discretization remains
numerical. PASS only covers the named checks, not the analytical theorem.

MANIFEST.json hashes package files except itself. Compare those hashes
before replay; replay may change floating-point output on another platform.
