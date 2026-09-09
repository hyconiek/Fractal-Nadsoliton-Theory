# Thirty-round FIN campaign: audit of a proposed learning rule

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
11. ST8631 — exact disjoint-minor rank bounds for strict and both declared
    legacy carriers; missing diagonals cannot evade the covariance bound.
12. ST8632 — complete stationary strict density blocks and exact rank bounds:
    seven in the interior gamma range, six at its PSD endpoint.
13. ST8633 — rank/entropy/purity extrema retain 32 phase choices; conservation
    of the actual state's spectrum prevents interpreting that as a selection law.
14. ST8634 — exact global energy-gap identity and attainable minima on each
    state-spectrum orbit.
15. ST8635 — a 50-dimensional local family of equally minimizing kernels.
16. ST8636 — a 39-dimensional family remains after fixing row sums and spectrum;
    a numerical curve checks, but does not substitute for, the local theorem.
17. ST8637 — exact circulant census leaves only two positive assignments,
    one C12 automorphism orbit; stronger premises genuinely narrow the answer.
18. ST8638 — a rigorous trajectory bound proves that the raw update can leave
    the positive strict-rate cone.
19. ST8639 — tangent-cone projection repairs viability but changes the law;
    a sparse pure equilibrium does not close the strict source problem.
20. ST8640 — dynamic shared spectral structure also requires regularity;
    the actual update can lose it while all weights are still positive.
21. ST8641 — exact self-reproduction from an encoded mixed-state seed;
    distinguish existence of a trajectory from derivation of its input.
22. ST8642 — matched Dirichlet-gradient propagation/learning rules; positivity
    and a conditional-negative-definiteness obstruction for strict.
23. ST8643 — viable regular positive projected learning, plus an exact pure
    spatial-averaging escape; the changed source map remains an explicit premise.
24. ST8644 — precision-domain audit and an exact positive weighted-graph
    Green parent after restoring a diagonal, with remaining parent freedom.
25. ST8645 — precision-flow sign audit, finite boundary failure of ascent and
    globally convergent fixed-covariance descent.
26. ST8646 — full chain rule and a legitimate joint Gaussian bootstrap;
    reciprocal consistency gives a flat operator family, not selection.
27. ST8647 — scalar Gibbs bootstrap has at most two spectral roots and cannot
    source seven strict sectors; small-decay near-zero residuals are not proof.
28. ST8648 — fast-learning tracking bound and Hamiltonian reduced limit;
    rapid tracking need not be rapid selection.
29. ST8649 — normal-ordering covariance fixes a diagonal algebra; abstract
    isospectrality is not automatically gauge equivalence of the learning model.
30. ST8650 — adversarial synthesis and requirement-by-requirement completion
    audit; no claim of completing FIN as a physical theory.

The thirty rounds are documented in REPORT.md. The final request audit and
test/replay verification must succeed before the goal is marked complete.
The previous campaign and existing archival sources remain unchanged;
completion of this research campaign does not mean physical closure of FIN.

## Replay

```sh
python3 fin_projected_learning/research.py
python3 fin_projected_learning/geometry.py
python3 fin_projected_learning/completion.py
python3 -m unittest discover -s fin_projected_learning -p 'test_*.py' -v
python3 fin_projected_learning/verify.py
```

There are 51 scientific tests, a controlled 30,000-step replay of the
original implementation, and exact rational separation of the seven strict
Fourier-sector eigenvalues. The existing weight-enclosure checker is reused
from `fin_replication_consistency/certify.py`; no old result is counted as a
new research round merely because it was replayed.
The new determinant and 120-assignment certificates use exact outward
arithmetic. Manifold dimensions rely on the accompanying rank proofs, not
only floating singular-value counts.
