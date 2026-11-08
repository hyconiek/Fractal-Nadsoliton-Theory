# Replay

From the repository root with NumPy and SymPy:

    OPENBLAS_NUM_THREADS=1 python fin_pde_two_mode_20261108/two_mode.py
    python fin_pde_two_mode_20261108/analyze.py

The solver takes several minutes and writes all half-unit diagnostics.
The analysis reconstructs fixed theoretical coefficients, two fit windows,
and reports both mixed-action estimates. Auxiliary initial preparations
are stated in REPORT_EN.md and do not replace the historical AW data.
The canonical.py dependency is imported from the frozen pilot package.
