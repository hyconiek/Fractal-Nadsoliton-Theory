# Reproduce this proof-oriented handoff

Python >=3.10 with numpy, scipy and mpmath installed (for diagnostics). The exact-rational certificate additionally uses only Python's standard library.

From this directory run:

```bash
python check_fold_chart.py
python SF-004/first_passage_certificate.py
python SF-004/independent_first_passage_check.py
python independent_verify.py
python verify_package.py
```

`check_fold_chart.py` regenerates binary64 `SF-001/NUMERICAL_CHART.json` and `first_passage_certificate.py` regenerates the exact-rational JSON, which are inputs to the independent verifier. On different numpy/scipy environments nominal numerical digits may drift; certificate status is algebraic/fraction arithmetic only. To validate archived hashes, run `sha256sum -c MANIFEST.sha256` **before** running writers, or on a clean extracted copy. `verify_package.py` checks file completeness and byte hashes and does not transform scientific outputs. The self-contained human proof is the set of Markdown documents; scripts verify arithmetic identities/inputs, not functional-analytic tightness or stochastic stopping theorems.

The earlier TH handoff is preserved in `source_context/` unchanged. Source R7P-031 exact-fold certificate is referenced by its public pinned URL and prior intake; its raw original bytes were not materialized in this execution, so its current hash cannot be claimed from this environment.
