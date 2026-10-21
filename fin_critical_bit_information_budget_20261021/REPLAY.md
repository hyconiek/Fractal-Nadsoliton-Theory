# Exact bounded replay

From the repository root:

```sh
prlimit --as=1073741824 --cpu=55 -- timeout 55s env PYTHONDONTWRITEBYTECODE=1 python3 fin_critical_bit_information_budget_20261021/verify_research.py --check --output /tmp/fin-bi-20261021.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_strategy_audit_20261019/verify_strategy.py --inputs-only
git diff --check
```

The new verifier uses standard-library Fraction arithmetic and rational
logarithm enclosures, not binary64 entropy equalities. It derives the actual
q12/N2 g0 collective-parity-count skeleton at model lag log2, checks a
genuine non-Markov binary collapse, finite-block lower/Fano/continuity
controls, analytical constant inequalities and a small de Bruijn source.
These are diagnostics, not a critical-domain theorem at N2, a native
small-lag uniform error or a mechanical proof of the boundary-layer limit.

Inputs and old listed manifest members are checked before/after; historical
writers and source finalizers are not run. --record authors only this fresh
package's metadata; --check writes only its optional /tmp output. No
external package, browser, PDF cache or experiment is a replay dependency.
