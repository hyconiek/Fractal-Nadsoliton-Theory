# Replay

From the repository root with Python and SymPy:

    python fin_core_state_map_bridge_20261113/check_bridge.py
    python fin_core_state_map_bridge_20261113/verify.py

The first command overwrites only this package's deterministic results.json.
The second verifies hashes and repeats the symbolic computation, requiring
byte-identical results. It does not run historical QW generators or rewrite
historical scientific packages. No simulation, external dataset, fitted
coefficient, or Lean toolchain is needed. Fetch any absent historical inputs
at the exact parent commit recorded in INPUTS.json; this workspace originally
contained only a partial repository mirror.

The proof in REPORT_EN.md treats a constant linear state map, common canonical
kinetics and local polynomial forces. It does not cover arbitrary nonlinear
reductions, effective integration of modes, or stability of the embedding.
