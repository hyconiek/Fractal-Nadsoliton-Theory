# Bounded replay and package check

From repository root, with sympy, numpy and scipy already installed:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 timeout 60s python3 fin_rp004_flux_review_20261019/verify_intake.py --replay --output /tmp/fin-rp004-review-check.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_rp004_flux_review_20261019/finalize_review.py --check
PYTHONDONTWRITEBYTECODE=1 python3 fin_intrinsic_entropy_source_review_20261019/finalize_review.py --check
PYTHONDONTWRITEBYTECODE=1 python3 fin_strategy_audit_20261019/verify_strategy.py --inputs-only
```

The first command copies source to /tmp before running both inspected
scientific scripts. Each child has a 60-second CPU/wall and 1-GiB address-
space bound. The outer timeout bounds this new verifier as well. No PDE
simulation or heavy algebra is run. Ordinary spectral quadrature may differ
in final digits across library environments; no byte identity is promised.
Its eight fixtures use momentum cutoff 18 and pay no outward integral-tail
bound. The exact symbolic/Fraction controls are recorded by name.

The root source `sha256sum -c MANIFEST.sha256` **does not pass**: finalize.py
is absent. The two task manifests do pass. The new verifier intentionally
reports that known omission separately from any hash mismatch. Do not run
or recreate an absent historical finalizer to manufacture a complete source.

The new finalizer checks metadata/hashes/links/JSON/AST and preservation,
not the truth of a PDE theorem. With no flag it seals **this new review only**;
`--check` is read-only. It does not edit source, old manifests, old AGENTS
hashes, the 47-file planning snapshot or the prior ID overlay. No internet
reading or imported literature theorem is used by the new elementary
pole/stationary-phase proof. Existing Darboux/NS inputs are explicitly reused.
