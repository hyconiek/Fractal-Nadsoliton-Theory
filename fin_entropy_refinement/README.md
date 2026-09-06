# FIN ST8547 — Entropy, binary refinement and hidden correlation

This package contains one analytical study with explicit proofs and 15
numerical tests. It follows ST8397–ST8546.

Main findings:

- Binary mobility compatibility alone permits a positive log-periodic entropy
  family and does not select Shannon.
- Binary compatibility plus diagonal concavity of mobility (implied by convex
  flux action) selects Shannon and logarithmic-mean mobility, up to entropy
  normalization, in the stated edge-local mass-coordinate class.
- Binary and ternary compatibility is an alternative sufficient premise.
- Coarse Onsager and entropy-production deficits are exactly controlled by
  correlations between base and fiber; equality holds for product states on
  a connected positive graph.

Read report.pdf or PROOFS_EN.md for complete hypotheses, proofs and limits.
The entropy selection is conditional. It does not derive FIN's kernel,
physical units, laboratory apparatus or a Theory of Everything.
The logarithmic-mean heat-gradient representation was already established
by Jan Maas in 2011.

Run from this directory:

    python3 analysis.py --output results.json
    python3 -m unittest -v test_analysis.py
    MPLCONFIGDIR=/tmp/fin-st8547-mpl python3 plot_results.py
    pdflatex -interaction=nonstopmode -halt-on-error report.tex
    pdflatex -interaction=nonstopmode -halt-on-error report.tex

The scientific tests recompute mathematical quantities. The checksum manifest
only protects artifact integrity. Neither constitutes laboratory evidence.

Creator: Krzysztof Żuchowski  
ORCID: 0009-0002-0909-3613  
Affiliation: Independent Researcher — Fractal Information Theory Project  
Resource type: Publication — Preprint  
Language: English  
Date: 2026-09-05  
Version: ST8547-v1  
License: CC BY 4.0 for the report; no new DOI is claimed or assigned.
