# KR intake replay

From repository root with numpy/scipy available for the versioned native
diagnostic; finite rational checks and producer replay use the stdlib:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 fin_stationary_kinetic_rigidity_review_20261016/verify_intake.py --replay --output /tmp/fin-kr-intake-replay.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_stationary_kinetic_rigidity_review_20261016/finalize_review.py --check
```

The intake verifier does not import the incoming producer. It independently
constructs a different nonuniform rational Gibbs chain, checks native local
characteristics, and copies the source before running its package checker
and scientific writer. Old archives, proofs and source JSON remain untouched.
verification.json stores actual logs and byte status. Its /tmp directory
is replay metadata, not a persistent dependency.

The original source is46 files plus its manifest. Five of six original
references are matched locally; old AGENTS and historical access remain
unverified. Actual current review inputs and older AGENTS prefixes are
bound separately. Do not reconstruct missing remote custody from hashes.

The package checker is mechanical. --check is read-only; without it only
this review's PACKAGE_CHECKS.json and MANIFEST.sha256 are generated.
