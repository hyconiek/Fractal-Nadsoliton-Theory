# UF intake replay and independent proof obligations

From repo root:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 fin_uniform_critical_review_20261014/check_gibbs_certificate.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 fin_uniform_critical_review_20261014/verify_intake.py --replay --output /tmp/fin-uf-intake.json
```

Numpy/scipy are needed for diagnostics. The Gibbs certificate checker uses
only the standard library and actual hash-bound rational spectral inputs.
It checks all payments and rejects a deliberately stronger exponent.

The intake hashes61 original source entries,4 source and30 stage inputs,
and73 KI archive files. Scientific writers execute only on a fresh `/tmp` copy
without historical caches; source files remain unchanged. The source validator
checks repeatability of the regenerated numerical instance in its environment,
not cross-environment byte identity with the original producer. `make_reports`
is not executed.

Manually inspect UF-001–004 and the source-bound SF joint-clock completion
before reusing the accepted process theorem. UI/absorption and finite-horizon
localization are analytical steps, not automatically checked by scripts.
Source spectral uncertainty is inherited; neither a point gain nor diagnostic
eigenvalue is a replacement for the exact operator-defined critical sequence.

The MT proof route is only a research strategy. No test here licenses a local
mixing gap or an exponentially long finite-N effective-law window.
