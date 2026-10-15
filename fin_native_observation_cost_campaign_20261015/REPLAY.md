# NOC-only replay

1. `python verify_noc.py` recomputes independent finite checks into `INDEPENDENT_VERIFICATION.json` (57/57).
2. `python verify_package.py` checks the manifest against current files and prints the number of covered entries.
3. Read `NOC-001/UNIFORM_MOMENTS_AND_CHARACTERISTICS.md`, `NOC-002/UNIFORM_DERIVATIVE_BOUNDS.md`, `NOC-003/NATIVE_MARTINGALE_DECOMPOSITION.md`, `NOC-004/FINITE_N_CONSISTENCY_OR_BLOCKER.md`. Their stochastic-limit dependencies are not mechanically verified by the scripts.
4. Compare source file hashes only for files listed in `SOURCE_INPUTS.sha256`; no unavailable review proof input is pretended present.
5. Repack source by `python finalize_noc.py` on the same input paths, then `python verify_package.py`.

Important: this package does not include any historical review, archive, source code or AGENTS change beyond hash references.
