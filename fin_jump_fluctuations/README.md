# FIN ST8548 — Markov fluctuations and quadratic entropy transport

Read report.pdf for the publication manuscript and PROOFS_EN.md for the
full analytical derivations, assumptions, comparisons and corrections.

Main results:

- The same FIN heat equation has a quadratic logarithmic-mean representation
  and a nonquadratic cosh representation from independent Markov walkers.
- Their local variances differ by z coth(z), where
  z is half the log population ratio. The fourth cumulant distinguishes them
  even at equilibrium.
- For sufficiently unequal populations, the quadratic noise prescription
  violates the local variance bound for integer-valued Markov counts.
- The exact current action and covariance commute with product-fiber
  projection for correlated states; the quadratic Onsager metric generally
  does not.

These are conditional finite-model statements. The microscopic counting
model and its physical realization are not derived from FIN.
The basic generalized-gradient construction is known literature and is
attributed in the report.

Reproduce in this directory:

    python3 analysis.py --output replay_results.json
    python3 -m unittest -v test_analysis.py
    MPLCONFIGDIR=/tmp/fin-st8548-mpl python3 plot_results.py
    pdflatex -interaction=nonstopmode -halt-on-error report.tex
    pdflatex -interaction=nonstopmode -halt-on-error report.tex

The tests recompute the mathematics, including finite-time tilted-generator
moments. They do not substitute file hashes for scientific validation.
No network or large computation is needed for replay.

Creator: Krzysztof Żuchowski  
ORCID: 0009-0002-0909-3613  
Affiliation: Independent Researcher — Fractal Information Theory Project  
Publication — Preprint; language English; version ST8548-v1  
Date: 2026-09-06  
License: CC BY 4.0 for the report. No DOI has been assigned in this package.
