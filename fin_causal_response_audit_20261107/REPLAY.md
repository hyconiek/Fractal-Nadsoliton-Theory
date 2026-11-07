# Replay

From the repository root with NumPy and SymPy:

    python fin_causal_response_audit_20261107/check_response.py

Recorded versions: NumPy 2.3.5, SymPy 1.14.0. The script writes
verification.json after checking pinned inputs and all controls.
The base grid has 16384 points on [-512,512); refinements use 32768
points on that interval and 24576 on [-768,768). Time propagation uses
exact Fourier sine/cosine multipliers. Space is still discretized and
periodic. Reported near-roundoff convergence does not certify the continuum
bound; the analytical proof supplies its uniform statement.

MANIFEST.json covers all other package files. Historical files are read
only. Replay on another numerical platform may change final digits.
