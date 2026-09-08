# New thirty-round FIN goal: audit of a sourced learning rule

Started 8 September 2026 after the previous ST8591–ST8620 campaign was
completed. The goal was started again; no old round is counted twice.
Outputs remain Markdown and executable mathematical checks, with no PDFs.

The source is the actual `nadsoliton_neural_analysis.py` update, compared
with the kernel report's adaptive-law discussion and the carefully scoped
unprojected no-go in ST2208/ST2209. This is not a generic repeated bridge or
selector search.

## Current sequence (proofs in REPORT.md)

1. ST8621 — identify the implemented projected leaky-Hebb law, separating it
   from genuine Oja/PCA and unprojected covariance matching.
2. ST8622 — exact mean/source equation for the actual random-phase drive.
3. ST8623 — fixed-step fluctuations and a same-noise contraction test.
4. ST8624 — the twelve-state entropy target is mathematically unattainable.
5. ST8625 — preserve degenerate spectral blocks in self-generated covariance.
6. ST8626 — exact mixed-state projected self-consistent kernel family.
7. ST8627 — a pure-state time-average witness is not instantaneous stationarity.
8. ST8628 — Lyapunov identity and isospectral state invariants for the explicit
   coupled law.
9. ST8629 — zero learning can coexist with persistent state motion.
10. ST8630 — dense-positive strict kernel: invariant-set rigidity and the
    resulting pure-state stationary obstruction.

Rounds 11–30 are not complete. The next move must be determined by these
results and their falsification tests, not by filling a predetermined list.
The previous campaign and existing archival source files remain unchanged.

## Replay

```sh
python3 fin_projected_learning/research.py
python3 -m unittest discover -s fin_projected_learning -p 'test_*.py' -v
python3 fin_projected_learning/verify.py
```

There are 17 scientific tests, a controlled 30,000-step replay of the
original implementation, and exact rational separation of the seven strict
Fourier-sector eigenvalues. The existing weight-enclosure checker is reused
from `fin_replication_consistency/certify.py`; no old result is counted as a
new research round merely because it was replayed.
