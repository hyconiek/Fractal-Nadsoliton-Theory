# Bounded read-only replay

Date: 2026-10-21. From the repository root:

```sh
prlimit --as=1073741824 --cpu=55 -- timeout 55s env PYTHONDONTWRITEBYTECODE=1 python3 fin_threshold_defect_virial_research_20261021/verify_research.py --check
prlimit --as=1073741824 --cpu=55 -- timeout 55s env PYTHONDONTWRITEBYTECODE=1 python3 fin_quadratic_localization_research_20261021/verify_research.py --check
prlimit --as=1073741824 --cpu=55 -- timeout 55s env PYTHONDONTWRITEBYTECODE=1 python3 fin_two_mode_dissipation_research_20261021/verify_research.py --check
git diff --check
```

The standard-library Fraction script checks exact rational-function weight
identities, sum-of-squares bounds, the interval-mean Fubini kernel, BOTH
forced virial identities as differential-polynomial flux identities,
Darboux threshold normalization and every reported rational constant.
It uses neither floating eigenvalues nor wave simulation as proof.

The own CONTINUUM_COERCIVITY_CERTIFICATE concerns a continuum block WITH
the rank defect. It is not CD's frozen reference-work certificate and
neither is a full nonlinear/modal FGR certificate. The trace history
budget is not tested into existence by finite flags.

Inputs, protected source manifests, own coverage, Python/JSON syntax and
Markdown links/whitespace are checked separately from scientific algebra.
`--record` authors only this fresh packet's metadata. Run it to completion
before a read-only check; never reseal a predecessor. No new browser,
PDF/cache dependency, independent agent, heavy algebra or external action.
