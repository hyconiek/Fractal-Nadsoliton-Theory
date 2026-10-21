# Read-only replay

Date: 2026-10-21. From the repository root:

```sh
prlimit --as=1073741824 --cpu=55 -- timeout 55s env PYTHONDONTWRITEBYTECODE=1 python3 fin_even_continuum_threshold_research_20261021/verify_research.py --check
git diff --check
```

The standard-library script uses Fraction polynomials. It checks even/odd
Jost numerator reductions, exact threshold moments, the null combination's
Taylor coefficient, actual HF and mixed moments, carrier-gap rational
bounds and negative/nonzero controls. pi and radicals are kept formal;
there is no floating near-zero, quadrature proof or large state enumeration.
The source identities are exact, while the analytical tail/modulation
proofs require their arguments in PROOF_OR_OPEN_LEMMA.md.

RESULTS reports the actual flag count separately from protected-manifest
entries and inputs. The script checks own manifest coverage, JSON/Python
syntax, Markdown links and whitespace. Hash/replay checks establish
preservation/reproducibility, not peer review or historical source custody.

`--record` authors only this fresh packet's metadata. Never reseal a
predecessor or rewrite its results to make a check pass. The earlier CM
packet is independently rechecked with its own read-only verifier.
No browser/network/cache, external agent, heavy algebra or physical action.
