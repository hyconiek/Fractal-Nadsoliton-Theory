# Replay

```sh
PYTHONDONTWRITEBYTECODE=1 python3 fin_balanced_phase_current_review_20261017/finalize_review.py --check
PYTHONDONTWRITEBYTECODE=1 python3 fin_balanced_phase_current_review_20261017/verify_intake.py --replay --output /tmp/fin-pc-intake-replay.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python3 fin_q12_physical_bridge_research_20261017/verify_bridges.py --output /tmp/fin-q12-physics-replay.json
```

Original source writers run only in a fresh /tmp copy. A derived wrapper
tests the original ZIP checker without claiming original-wrapper custody.
Independent full manifest coverage also binds the excluded checker and
metadata. build_handoff.py is not executed.

The new verifier uses exact rational controls and a complete q12/N3 iid
stratum-coupling check, not critical MT enumeration. Current predecessor
scientific inputs and historical AGENTS prefixes are checked separately.
The physical research has its own symbolic/finite/nominal scopes and
historical browser-capture log. Replaying math does not repeat the browser
history or validate original PDF hashes.
