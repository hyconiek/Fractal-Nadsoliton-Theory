# Replay

Python with NumPy, SymPy and Matplotlib is required. Recorded NumPy 2.3.5
and SymPy 1.14.0. Run from the repository root:

    python fin_canonical_audit_pde_pilot_20261107/check_canonical.py
    python fin_canonical_audit_pde_pilot_20261107/check_bridge.py
    OPENBLAS_NUM_THREADS=1 python fin_canonical_audit_pde_pilot_20261107/pde_pilot.py
    python fin_canonical_audit_pde_pilot_20261107/summarize_pilot.py

The PDE pilot takes several minutes. pde_results.json is a generated local
checkpoint containing all diagnostics and is excluded from version control.
The committed pilot_summary.json contains all run configurations and metrics;
trajectories.csv retains every fourth recorded diagnostic and each final
point. The solver and summarizer reproduce the full checkpoint and outputs.
The figure uses the full recorded diagnostics, not only the decimated CSV.

The code reconstructs G4 symbolically once per process. The inverse auxiliary
flows use 16 RK4 steps each; check_canonical compares 16 with 32 steps.
The time step is adjusted downward to end exactly at T. No result is labelled
an independent-team replication or empirical measurement. Manifest hashes
cover committed package files excluding the manifest itself.
