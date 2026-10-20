# Replay

From repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 fin_post_rp004_history_review_20261020/verify_intake.py --replay --output /tmp/fin-post-history-replay-20261020.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_post_rp004_history_review_20261020/finalize_review.py --check
PYTHONDONTWRITEBYTECODE=1 python3 fin_rp004_flux_review_20261019/finalize_review.py --check
PYTHONDONTWRITEBYTECODE=1 python3 fin_intrinsic_entropy_source_review_20261019/finalize_review.py --check
PYTHONDONTWRITEBYTECODE=1 python3 fin_strategy_audit_20261019/verify_strategy.py --inputs-only
```

The new verifier checks source hashes first, copies the complete bundle to
/tmp, and runs the master and five inspected scientific scripts there.
Each child has a 60-second CPU/wall and 1-GiB address-space ceiling. The KL
script's one hardcoded output path is mechanically relocated only in that
copy. Scientific expressions are not changed. A resource stop, if one
occurs in another environment, is not a mathematical counterexample.

Four numeric outputs differ here from incoming bytes; max numerical drift
is about 3.82e-14. Only the anonymous-history JSON is byte-identical in this
replay. Library-version changes are recorded separately. Exact Q(sqrt(3))
controls in this review use a different implementation, not producer helpers.
Other controls are ordinary floating numerical diagnostics, not intervals.

The full consolidated manifest covers 75 entries. Its own file and verifier
are independently bound; no absent ZIP is reconstructed or CRC-verified.
Original INPUTS paths are resolved by explicit mapping to local/extracted
counterparts with archived-only matches and missing originals distinguished.

The finalizer writes only this new intake's binding and manifest when run
without `--check`; with `--check` it is read-only. It does not rewrite sources,
AGENTS or earlier scientific/status packages. A pass certifies these package
checks, not the truth of every analytical theorem or a physical realization.
