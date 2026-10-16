# Source-preserving FP intake replay

From repository root, with mpmath available for separate numerical diagnostics:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 fin_count_partition_feasibility_review_20261016/verify_intake.py --replay --output /tmp/fin-fp-intake-replay.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_count_partition_feasibility_review_20261016/finalize_review.py --check
```

The intake checker imports no incoming scientific code. It independently
checks exact primal/dual, Fourier direction and rational threshold/cap bounds;
its erfc formula is a separate noncertified numerical preview.

With --replay it copies the source to /tmp, runs verify_package.py before
verify_fp.py, and records actual output byte status. Producer files are
checked again afterwards. Four archived source members are checked inside
the actual MT/UF ZIPs; missing BA/KR original wrappers are not reconstructed.

Current BA/KR/AC/MT/UF inputs and saved science hashes are independently
bound. Historical AGENTS snapshots are not equated with the current ledger;
old append-only input hashes use their actual matching prefixes where present.

The finalizer is mechanical and read-only with --check. Without that flag,
only this review's PACKAGE_CHECKS.json and MANIFEST.sha256 are generated.
It never edits old proofs, reports, producers or empirical evidence.
