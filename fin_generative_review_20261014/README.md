# FIN generative campaign: repository intake and architectural checkpoint

Review date: 2026-10-14. Baseline: `9e507988`.

[INTAKE_20261014.md](INTAKE_20261014.md) is the authoritative scope of
acceptance of [the supplied campaign](../fin_generative_campaign_20261014/HANDOFF.md).
[RG_C06_DECISION.md](RG_C06_DECISION.md) selects the next research branch.
[CLAIM_REGISTER.json](CLAIM_REGISTER.json) distinguishes exact identities,
validated interval results, numerical results and unfinished obligations.

Run the repository intake without modifying source artifacts:

```bash
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 fin_generative_review_20261014/verify_intake.py
```

The checker verifies the 63-entry source manifest, binds the reconstruction
to the versioned safe-seed model, independently encloses all 1,728 static
log-weights, computes rational ANOVA and rebuilds the 1,728-state generator.
Dynamic probabilities are checked with a Poisson uniformization series,
separately from the producer's Krylov propagation. Dynamic arithmetic is
still numerical, not an interval certificate.

`verification.json` records this intake run. Its interval endpoints are
serialized outward; the public conservative static bound is
`TV <= 0.051301995505` for the specified exact-decimal FIN model.

The supplied producer, interval checker and 14-check second verifier were
also run successfully in `/tmp/fin-gen-intake-YXUvh2ID/campaign`. Their
regenerated runtime values and historical output documents were not copied
over the immutable source package. `finalize_package.py` was not needed and
was not run: it depends on a historical uploaded-plan path.

Dependencies: Python 3, NumPy, SciPy and mpmath. No network, hardware,
FIN-M0 trajectory or unopened task335 gain is required by this checker.

Next execution plan:
[FIN_GENERATIVE_MEMORY_NEXT_CHAIN_20261014.md](../FIN_GENERATIVE_MEMORY_NEXT_CHAIN_20261014.md).
