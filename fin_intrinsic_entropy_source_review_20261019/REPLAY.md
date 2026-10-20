# Source-preserving entropy intake replay

From the actual repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 fin_intrinsic_entropy_source_review_20261019/verify_intake.py --replay --output /tmp/fin-id001-intake-replay.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_intrinsic_entropy_source_review_20261019/finalize_review.py --check
PYTHONDONTWRITEBYTECODE=1 python3 fin_strategy_audit_20261019/verify_strategy.py --inputs-only
git diff --check
```

Do not run the producer writer in its archived source directory. The intake
copies it to /tmp and reports floating drift and source hashes before/after.
The fixed native p(g3), quotient and entropy controls are point numerics,
not a finite-H provider, numerical KS computation or outward kernel proof.
Analytical entropy/decimation proofs are separate from finite controls.

Use only --inputs-only on the live historical strategy checker. Its full
mode authors VERIFICATION.json for the original UNEXECUTED architecture
snapshot; it is not a live campaign-status recorder. Do not rehash its
47-file manifest or rewrite the old task READY fields. Current execution
status is the separate FIN_RESEARCH_EXECUTION_STATUS_20261019.json overlay.

Historical AGENTS inputs are checked on their actual hash-bound prefixes;
the new current context is bound separately. Read the live instructions,
not a historical replay file. Finalizer write mode authors only this new
review's index/manifest. No experiment, source acquisition or Git write.
