# Replay

Use Python with NumPy and SymPy, from the repository root:

    OPENBLAS_NUM_THREADS=1 python fin_mixed_third_amplitude_20261113/run_case.py 0
    OPENBLAS_NUM_THREADS=1 python fin_mixed_third_amplitude_20261113/run_case.py 1
    OPENBLAS_NUM_THREADS=1 python fin_mixed_third_amplitude_20261113/run_case.py 2
    python fin_mixed_third_amplitude_20261113/analyze.py
    python fin_mixed_third_amplitude_20261113/verify.py

Optional figure regeneration requires Matplotlib:

    python fin_mixed_third_amplitude_20261113/plot.py

Cases are independent and can run concurrently. Existing phase*.json
checkpoints are skipped. Recompute in a fresh sibling directory containing
run_case.py or explicitly remove only the corresponding checkpoint first.
CSV files retain every (t,I,J) diagnostic. JSON stores run parameters and
maximum energy/flux residuals; spatial fields and full energy traces are not
archived. The previous campaign.py solver and compare.py estimator are
imported unchanged, with hashes pinned in INPUTS.json. Their frozen upstream
canonical.py and CD coefficient implementation are also pinned.

The refinement changes dx and dt together. It tests their combined effect,
not separate spatial and temporal convergence orders. No extrapolation of
solver error or rigorous continuum-limit bound is claimed. Fits compare the
same physical diagnostic nodes and slow-time windows as the previous study.
