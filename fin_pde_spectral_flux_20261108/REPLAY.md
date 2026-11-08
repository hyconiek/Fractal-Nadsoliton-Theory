# Replay

From the repository root, using NumPy and SymPy:

    OPENBLAS_NUM_THREADS=1 python fin_pde_spectral_flux_20261108/spectral_flux.py
    python fin_pde_spectral_flux_20261108/analyze.py

Six runs take several minutes. Results include diagnostic samples every five
model-time units, plus the final point; summary residuals use all half-unit
diagnostics. The canonical map is imported from the frozen pilot package.
The previous campaign is read only for the comparison. All values use
dimensionless model units. Checks are numerical quality gates, not proofs.
