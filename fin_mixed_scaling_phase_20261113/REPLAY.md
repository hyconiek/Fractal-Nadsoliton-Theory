# Replay

From the repository root, with Python, NumPy and SymPy:

    OPENBLAS_NUM_THREADS=1 python fin_mixed_scaling_phase_20261113/campaign.py
    python fin_mixed_scaling_phase_20261113/compare.py
    python fin_mixed_scaling_phase_20261113/verify.py

The campaign skips completed per-run JSON checkpoints. To recompute, copy
campaign.py into a fresh sibling directory (so the canonical import still
resolves), then compare the resulting CSV and metadata within numerical
roundoff. To regenerate in place, explicitly remove only the four eps*.json
checkpoints first. Each CSV retains every dense (t,I,J) diagnostic. Metadata
retains maxima of global energy drift and two physical local flux residuals;
full spatial fields and individual energy samples are not retained.

No fitted coefficient is fed back into the PDE. Analysis uses the frozen CD
self-channel coefficients and a free intercept. Integration-grid decimation
uses identical regression nodes and includes the final endpoint. The first
sample at or after each slow-time window boundary is used; its actual time
may differ from the requested boundary by less than one diagnostic interval.

verify.py checks reproducible numerical integrity and input/output hashes;
it does not prove the PDE action law or provide independent peer review.
