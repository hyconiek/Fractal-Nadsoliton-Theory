# FIN — Replication Consistency and Collective Activity

Programmes ST8579–ST8590. English report and local replay package.

## Principal results

For a strongly projectively consistent, exchangeable family of finite Markov
generators, define b_ij as the rate at which both coordinates change from
initial vertices (i,j). Unlimited population replication forces B=(b_ij) to be
a completely positive real matrix (a Gram matrix of nonnegative vectors).
This is not a statement about quantum completely positive maps.

Consequently b_ii=0 for every origin forces the independent tensor generator.
If all b_ii <= delta, the previous transition-TV theorem improves to
min(1, N(N-1) delta t). Unlimited replication and the diagonal bound remain
additional premises not derived from the existing FIN kernel.

An explicit strict model is consistent through three copies, reversible and
D12-symmetric, with zero same-origin activity but nonzero mixed activity. It
cannot extend consistently to four copies. Even maximizing all b_ij leaves
distinct irreducible models with the same singleton laws and equilibrium.

The neural audit encloses the actual transcendental weights with exact
rational intervals and classifies all 4096 threshold states. Certified counts:

| Declared synchronous network | Fixed points | Two-cycles |
|---|---:|---:|
| Strict cycle | 2 | 111 |
| Canonical legacy cycle | 6 | 64 |
| Canonical legacy line | 6 | 18 |

These are conditional spin-network results, not physical vacua, learned
memory capacities, or evidence for a ToE. A failed uniform-four-state test
and its analytic exception are explicitly retained.

## Reproduction

Dependencies: Python 3, NumPy, SciPy, mpmath and Matplotlib; pdflatex to build
the English PDF. Scripts install nothing and use no network.

From this directory:

```sh
python3 analysis.py
python3 certify.py
python3 -m unittest discover -s . -p 'test_*.py' -v
python3 build.py --verify --pdf --package
```

The final command records actual test output, checks exact replay equality,
renders figures, compiles the report twice and creates a SHA-256 manifest
and portable ZIP. A standalone replay preserves the supplied original
provenance manifest and writes a separate replay manifest.

The proof source is `report.tex`; `report.pdf` is the English manuscript.
`results.json` contains floating computations and the retained failed test.
`sign_certificate.json` contains exact outward integer weight tables and
the resulting certified successor-graph counts. `verification.json` records
executed test/replay outcomes. Hashes certify identity, not scientific truth.

Author: Krzysztof Żuchowski; ORCID 0009-0002-0909-3613.
Affiliation: Independent Researcher — Fractal Information Theory Project.
Report license: CC BY 4.0. No DOI assigned; no upload or laboratory claim.
