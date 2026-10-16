# MM intake replay

From the repository root, with numpy, scipy, sympy and mpmath installed:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 fin_mediator_memory_observability_review_20261016/certify_margin.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 fin_mediator_memory_observability_review_20261016/verify_intake.py --replay --output /tmp/fin-mm-intake-replay.json
PYTHONDONTWRITEBYTECODE=1 python3 fin_mediator_memory_observability_review_20261016/finalize_review.py --check
```

The intake verifier checks original hashes before and after, copies the
source to a temporary cache-free directory and executes verify_package.py
before verify_mm.py. It does not run build_mm_documents.py, reconstruct
missing uploaded wrappers or write producer files. verification.json stores
actual replay logs, source comparison and the temporary path; that path
need not persist in a clean checkout.

The certificate command prints a reconstructed object; it does not overwrite
the public certificate unless --output is explicitly supplied. Its backend
is mpmath.iv, with exact outward decimal conversion. The ordinary floating
matrix-semigroup diagnostic is separate and not an interval proof.

finalize_review.py --check is read-only. Without --check it authors only
this review's mechanical PACKAGE_CHECKS.json and MANIFEST.sha256 after
checking actual dependencies, source scope and the conditional plan.
